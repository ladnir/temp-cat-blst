use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    if env::var("CARGO_FEATURE_SERDE_SECRET").is_ok() {
        println!(
            "cargo:warning=blst: non-production feature serde-secret enabled"
        );
    }

    // account for cross-compilation [by examining environment variables]
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap();
    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap();
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let target_no_std = target_os.eq("none")
        || (target_os.eq("unknown") && target_arch.eq("wasm32"))
        || target_os.eq("uefi")
        || env::var("BLST_TEST_NO_STD").is_ok();

    if !target_no_std {
        println!("cargo:rustc-cfg=feature=\"std\"");
        if target_arch.eq("wasm32") || target_os.eq("unknown") {
            println!("cargo:rustc-cfg=feature=\"no-threads\"");
        }
    }
    println!("cargo:rerun-if-env-changed=BLST_TEST_NO_STD");

    /*
     * Use pre-built libblst.a if there is one. This is primarily
     * for trouble-shooting purposes. Idea is that libblst.a can be
     * compiled with flags independent from cargo defaults, e.g.
     * '../../build.sh -O1 ...'.
     */
    if Path::new("libblst.a").exists() {
        println!("cargo:rustc-link-search=.");
        println!("cargo:rustc-link-lib=blst");
        println!("cargo:rerun-if-changed=libblst.a");
        return;
    }

    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let blst_base_dir = resolve_blst_source_dir(&manifest_dir);
    println!("Using blst source directory {}", blst_base_dir.display());
    println!(
        "cargo:rerun-if-changed={}",
        blst_base_dir.join("CMakeLists.txt").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        blst_base_dir.join("blst").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        blst_base_dir.join("bindings").display()
    );

    let portable = cfg!(feature = "portable");
    let enable_avx = cfg!(feature = "avx");
    let enable_avx_portable = enable_avx && target_os.eq("windows");
    let enable_avx_native = enable_avx && !enable_avx_portable;
    if portable && cfg!(feature = "force-adx") {
        panic!("cannot compile with both `portable` and `force-adx` features");
    }
    if portable && enable_avx {
        panic!("cannot compile with both `portable` and `avx` features");
    }

    let enable_asm = !portable
        && !target_no_std
        && (target_arch.eq("x86_64") || target_arch.eq("aarch64"));
    let enable_adx = enable_asm && should_enable_adx(&target_arch, &target_env);
    let profile = env::var("PROFILE").unwrap_or_else(|_| "debug".to_string());
    let build_type = if profile.eq_ignore_ascii_case("release") {
        "Release"
    } else {
        "Debug"
    };
    let build_dir = PathBuf::from(env::var("OUT_DIR").unwrap()).join("cmake");
    let use_ninja =
        env::var_os("CMAKE_GENERATOR").is_none() && command_works("ninja", &["--version"]);

    let mut configure = Command::new("cmake");
    configure.arg("-S").arg(&blst_base_dir);
    configure.arg("-B").arg(&build_dir);
    if use_ninja {
        configure.args(["-G", "Ninja"]);
    }
    configure.arg(format!("-DCMAKE_BUILD_TYPE={build_type}"));
    if target_env.eq("msvc") {
        configure.arg("-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL");
        if let Some(cl_path) = find_msvc_cl() {
            configure.arg(format!("-DCMAKE_C_COMPILER={}", cl_path.display()));
            configure.arg(format!("-DCMAKE_CXX_COMPILER={}", cl_path.display()));
        }
    }
    configure.args([
        "-DBLST_BUILD_SHARED=OFF",
        "-DBLST_TEST=OFF",
        "-DBLST_BENCHMARK=OFF",
        "-DBLST_CUDA=OFF",
        "-DBLST_CUDA_MSM=OFF",
        "-DBLST_PARALLEL_MSM=OFF",
        "-DBLST_CRNS_BATCH_VROOM_VENDOR=OFF",
    ]);
    configure.arg(format!("-DBLST_PORTABLE={}", on_off(portable)));
    configure.arg(format!("-DBLST_ASM={}", on_off(enable_asm)));
    configure.arg(format!("-DBLST_DEFAULT_ASM={}", on_off(enable_asm)));
    configure.arg(format!("-DBLST_ADX={}", on_off(enable_adx)));
    configure.arg(format!("-DBLST_AVX512={}", on_off(enable_avx_native)));
    configure.arg(format!("-DBLST_AVX512_PORTABLE={}", on_off(enable_avx_portable)));
    configure.arg(format!("-DBLST_DEFAULT_AVX={}", on_off(enable_avx)));
    run(&mut configure, "configure blst with CMake");

    let mut build = Command::new("cmake");
    build.arg("--build").arg(&build_dir);
    build.args(["--target", "blst_static"]);
    if !use_ninja && target_env.eq("msvc") {
        build.args(["--config", build_type]);
    }
    run(&mut build, "build blst_static");

    let lib_name = if target_env.eq("msvc") {
        "blst.lib"
    } else {
        "libblst.a"
    };
    let lib_dir = find_parent_dir_containing(&build_dir, lib_name)
        .unwrap_or_else(|| panic!("unable to locate {lib_name} under {}", build_dir.display()));

    println!("cargo:rustc-link-search=native={}", lib_dir.display());
    println!("cargo:rustc-link-lib=static=blst");
    if enable_avx {
        emit_avx_link_deps(&blst_base_dir, &target_env, build_type, &build_dir);
    }
    if !target_env.eq("msvc") {
        let cpp_runtime = if target_os.eq("macos") || target_os.eq("ios") {
            "c++"
        } else {
            "stdc++"
        };
        println!("cargo:rustc-link-lib={cpp_runtime}");
    }

    // pass some DEP_BLST_* variables to dependents
    println!(
        "cargo:BINDINGS={}",
        blst_base_dir.join("bindings").to_string_lossy()
    );
    println!(
        "cargo:C_SRC={}",
        blst_base_dir.join("blst").to_string_lossy()
    );
}

fn on_off(value: bool) -> &'static str {
    if value { "ON" } else { "OFF" }
}

fn resolve_blst_source_dir(manifest_dir: &Path) -> PathBuf {
    let packaged_root = manifest_dir.join("blst");
    if packaged_root.join("CMakeLists.txt").is_file() {
        return packaged_root;
    }

    if manifest_dir.join("CMakeLists.txt").is_file() {
        return manifest_dir.to_path_buf();
    }

    manifest_dir
        .parent()
        .and_then(|dir| dir.parent())
        .filter(|dir| dir.join("CMakeLists.txt").is_file())
        .map(Path::to_path_buf)
        .unwrap_or_else(|| {
            panic!(
                "unable to locate blst source root from {}",
                manifest_dir.display()
            )
        })
}

fn should_enable_adx(target_arch: &str, target_env: &str) -> bool {
    if target_arch != "x86_64" {
        return false;
    }
    if cfg!(feature = "force-adx") || target_env.eq("sgx") {
        return true;
    }
    if env::var("CARGO_ENCODED_RUSTFLAGS")
        .unwrap_or_default()
        .contains("target-cpu=")
    {
        return env::var("CARGO_CFG_TARGET_FEATURE")
            .unwrap_or_default()
            .split(',')
            .any(|feature| feature == "adx");
    }
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("adx")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

fn command_works(command: &str, args: &[&str]) -> bool {
    Command::new(command)
        .args(args)
        .status()
        .map(|status| status.success())
        .unwrap_or(false)
}

fn find_msvc_cl() -> Option<PathBuf> {
    let from_where = Command::new("where.exe")
        .arg("cl.exe")
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| {
            String::from_utf8(output.stdout).ok().and_then(|stdout| {
                stdout
                    .lines()
                    .map(str::trim)
                    .find(|line| line.to_ascii_lowercase().ends_with("\\cl.exe"))
                    .map(PathBuf::from)
            })
        });
    if from_where.is_some() {
        return from_where;
    }

    let vswhere = PathBuf::from(r"C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe");
    if !vswhere.is_file() {
        return None;
    }

    Command::new(vswhere)
        .args([
            "-latest",
            "-products",
            "*",
            "-requires",
            "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
            "-find",
            r"VC\Tools\MSVC\**\bin\Hostx64\x64\cl.exe",
        ])
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| {
            String::from_utf8(output.stdout).ok().and_then(|stdout| {
                stdout
                    .lines()
                    .map(str::trim)
                    .find(|line| line.to_ascii_lowercase().ends_with("\\cl.exe"))
                    .map(PathBuf::from)
            })
        })
}

fn emit_avx_link_deps(
    blst_base_dir: &Path,
    target_env: &str,
    build_type: &str,
    cargo_build_dir: &Path,
) {
    if target_env.eq("msvc") {
        let vcpkg_root = blst_base_dir
            .parent()
            .expect("blst repo should have a parent directory")
            .join("vcpkg")
            .join("installed")
            .join("x64-windows");
        let lib_dir = if build_type.eq_ignore_ascii_case("release") {
            vcpkg_root.join("lib")
        } else {
            vcpkg_root.join("debug").join("lib")
        };
        let dll_path = if build_type.eq_ignore_ascii_case("release") {
            vcpkg_root.join("bin").join("mpir.dll")
        } else {
            vcpkg_root.join("debug").join("bin").join("mpir.dll")
        };
        println!("cargo:rustc-link-search=native={}", lib_dir.display());
        println!("cargo:rustc-link-lib=mpir");
        copy_runtime_dll(&dll_path, cargo_build_dir);
    } else {
        println!("cargo:rustc-link-lib=gmpxx");
        println!("cargo:rustc-link-lib=gmp");
    }
}

fn copy_runtime_dll(dll_path: &Path, cargo_build_dir: &Path) {
    println!("cargo:rerun-if-changed={}", dll_path.display());
    let Some(profile_dir) = cargo_build_dir.ancestors().nth(4) else {
        return;
    };
    let copies = [profile_dir.join("mpir.dll"), profile_dir.join("deps").join("mpir.dll")];
    for target in copies {
        if let Some(parent) = target.parent() {
            let _ = fs::create_dir_all(parent);
        }
        let _ = fs::copy(dll_path, target);
    }
}

fn run(command: &mut Command, context: &str) {
    let status = command
        .status()
        .unwrap_or_else(|err| panic!("failed to {context}: {err}", context = context, err = err));
    if !status.success() {
        panic!(
            "{context} failed with status {status}",
            context = context,
            status = status
        );
    }
}

fn find_parent_dir_containing(root: &Path, file_name: &str) -> Option<PathBuf> {
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let entries = fs::read_dir(&dir).ok()?;
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.file_name().and_then(|name| name.to_str()) == Some(file_name) {
                return path.parent().map(Path::to_path_buf);
            }
        }
    }
    None
}

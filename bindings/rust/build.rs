use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap();
    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap();
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let host = env::var("HOST").unwrap();
    let target = env::var("TARGET").unwrap();
    let cross_compiling = host != target;
    let target_no_std = target_os.eq("none")
        || (target_os.eq("unknown") && target_arch.eq("wasm32"))
        || target_os.eq("uefi")
        || env::var("BLST_TEST_NO_STD").is_ok();

    println!("cargo:rustc-check-cfg=cfg(blst_std)");
    if !target_no_std {
        println!("cargo:rustc-cfg=blst_std");
    }
    for name in [
        "BLST_TEST_NO_STD",
        "BLST_ALLOW_UNTESTED_ARM64",
        "BLST_ALLOW_UNTESTED_CHERI",
        "BLST_PREBUILT_LIB",
        "BLST_PREBUILT_CONFIG_INCLUDE",
        "CMAKE_GENERATOR",
        "CMAKE_TOOLCHAIN_FILE",
        "CC",
        "CXX",
        "CFLAGS",
        "CXXFLAGS",
        "BLST_FUZZ_NATIVE_SANITIZER",
        "CARGO_ENCODED_RUSTFLAGS",
    ] {
        println!("cargo:rerun-if-env-changed={name}");
    }

    let allow_untested_arm64 = env_flag("BLST_ALLOW_UNTESTED_ARM64");
    let allow_untested_cheri = env_flag("BLST_ALLOW_UNTESTED_CHERI");
    let target_is_cheri = env::var("CARGO_CFG_TARGET_POINTER_WIDTH").as_deref()
        == Ok("128")
        || target_env.to_ascii_lowercase().contains("cheri")
        || target_env.to_ascii_lowercase().contains("purecap");
    if target_is_cheri && !(allow_untested_arm64 && allow_untested_cheri) {
        panic!(
            "CHERI pure-capability builds are experimental and unqualified; set both \
             BLST_ALLOW_UNTESTED_ARM64=1 and BLST_ALLOW_UNTESTED_CHERI=1 to attempt one"
        );
    }
    if target_arch.eq("aarch64")
        && !target_os.eq("macos")
        && !allow_untested_arm64
    {
        panic!(
            "ARM64 builds are experimental and unqualified; set \
             BLST_ALLOW_UNTESTED_ARM64=1 to attempt one"
        );
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
    let force_adx = cfg!(feature = "force-adx");
    let opt_level = env::var("OPT_LEVEL").unwrap_or_else(|_| "0".to_string());
    let build_type = native_build_type(&opt_level);
    validate_build_policy(&target_arch, portable, enable_avx, force_adx)
        .unwrap_or_else(|message| panic!("{}", message));

    let enable_avx_portable = enable_avx && target_os.eq("windows");
    let enable_avx_native = enable_avx && !enable_avx_portable;

    let toolchain = env::var_os("CMAKE_TOOLCHAIN_FILE").map(PathBuf::from);
    validate_cross_toolchain(&host, &target, toolchain.as_deref())
        .unwrap_or_else(|message| panic!("{}", message));

    if let Some(prebuilt) = env::var_os("BLST_PREBUILT_LIB") {
        let prebuilt = PathBuf::from(prebuilt);
        let config_include = env::var_os("BLST_PREBUILT_CONFIG_INCLUDE")
            .map(PathBuf::from)
            .unwrap_or_else(|| {
                panic!(
                    "BLST_PREBUILT_CONFIG_INCLUDE is required with BLST_PREBUILT_LIB"
                )
            });
        link_prebuilt(
            &prebuilt,
            &config_include,
            &target_env,
            &target_os,
            &blst_base_dir,
        );
        return;
    }

    let enable_asm = !portable
        && !target_no_std
        && (target_arch.eq("x86_64") || target_arch.eq("aarch64"));
    let enable_adx = enable_asm
        && should_enable_adx(&target_arch, &target_env, cross_compiling);
    let user_generator = env::var("CMAKE_GENERATOR").ok();
    let use_ninja = user_generator.is_none()
        && !target_env.eq("msvc")
        && command_works("ninja", &["--version"]);
    let msvc_generator = if target_env.eq("msvc") && user_generator.is_none() {
        find_visual_studio_generator()
    } else {
        None
    };
    let build_dir_prefix = if use_ninja {
        "cmake-ninja"
    } else if target_env.eq("msvc") && user_generator.is_none() {
        "cmake-msvc"
    } else {
        "cmake"
    };
    let build_dir_name =
        format!("{}-{}", build_dir_prefix, build_type.to_ascii_lowercase());
    let build_dir =
        PathBuf::from(env::var("OUT_DIR").unwrap()).join(&build_dir_name);

    let mut configure = Command::new("cmake");
    configure.arg("-S").arg(&blst_base_dir);
    configure.arg("-B").arg(&build_dir);
    if use_ninja {
        configure.args(["-G", "Ninja"]);
    } else if target_env.eq("msvc") && user_generator.is_none() {
        if let Some(generator) = msvc_generator {
            configure.args(["-G", generator]);
        }
        if let Some(generator_arch) = msvc_cmake_arch(&target_arch) {
            configure.args(["-A", generator_arch]);
        }
    }
    if let Some(path) = toolchain.as_deref() {
        configure
            .arg(format!("-DCMAKE_TOOLCHAIN_FILE={}", path.to_string_lossy()));
    }
    configure.arg(format!("-DCMAKE_BUILD_TYPE={build_type}"));
    if target_env.eq("msvc") {
        configure.arg("-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL");
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
    configure.arg(format!(
        "-DBLST_ALLOW_UNTESTED_ARM64={}",
        on_off(allow_untested_arm64)
    ));
    configure.arg(format!(
        "-DBLST_ALLOW_UNTESTED_CHERI={}",
        on_off(allow_untested_cheri)
    ));
    configure.arg(format!("-DBLST_AVX512={}", on_off(enable_avx_native)));
    configure.arg(format!(
        "-DBLST_AVX512_PORTABLE={}",
        on_off(enable_avx_portable)
    ));
    configure.arg(format!("-DBLST_DEFAULT_AVX={}", on_off(enable_avx)));
    run(&mut configure, "configure blst with CMake");

    let mut build = Command::new("cmake");
    build.arg("--build").arg(&build_dir);
    build.args(["--target", "blst_static"]);
    if !use_ninja && target_env.eq("msvc") {
        build.args(["--config", build_type]);
    }
    run(&mut build, "build blst_static");

    let lib_dir = built_library_dir(&build_dir, &target_env, build_type);
    emit_native_link(&lib_dir, &target_env, &target_os);
    emit_fuzz_sanitizer_runtime(&target_env, &target_os);
    emit_metadata(&blst_base_dir, &build_dir.join("blst"));
}

fn native_build_type(opt_level: &str) -> &'static str {
    if opt_level == "0" {
        "Debug"
    } else {
        "Release"
    }
}

fn validate_build_policy(
    target_arch: &str,
    portable: bool,
    avx: bool,
    force_adx: bool,
) -> Result<(), &'static str> {
    if portable && force_adx {
        return Err(
            "cannot compile with both `portable` and `force-adx` features",
        );
    }
    if portable && avx {
        return Err("cannot compile with both `portable` and `avx` features");
    }
    if avx && target_arch != "x86_64" {
        return Err("the `avx` feature is supported only on x86_64 targets");
    }
    Ok(())
}

fn validate_cross_toolchain(
    host: &str,
    target: &str,
    toolchain: Option<&Path>,
) -> Result<(), String> {
    if host != target && toolchain.is_none() {
        return Err(format!(
            "cross-compiling CatBLST from {} to {} requires an explicit \
             CMAKE_TOOLCHAIN_FILE; refusing to compile a host archive for the target",
            host, target
        ));
    }
    if let Some(path) = toolchain {
        if !path.is_file() {
            return Err(format!(
                "CMAKE_TOOLCHAIN_FILE is not a file: {}",
                path.display()
            ));
        }
    }
    Ok(())
}

fn on_off(value: bool) -> &'static str {
    if value {
        "ON"
    } else {
        "OFF"
    }
}

fn env_flag(name: &str) -> bool {
    let Ok(value) = env::var(name) else {
        return false;
    };

    match value.to_ascii_lowercase().as_str() {
        "1" | "on" | "true" | "yes" => true,
        "" | "0" | "off" | "false" | "no" => false,
        _ => panic!(
            "{} must be one of: 1, ON, true, yes, 0, OFF, false, no",
            name
        ),
    }
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

fn should_enable_adx(
    target_arch: &str,
    target_env: &str,
    cross_compiling: bool,
) -> bool {
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
    if cross_compiling {
        return false;
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

fn find_visual_studio_generator() -> Option<&'static str> {
    let vswhere = PathBuf::from(
        r"C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe",
    );
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
            "-property",
            "installationVersion",
        ])
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| {
            String::from_utf8(output.stdout).ok().and_then(|stdout| {
                stdout.lines().map(str::trim).find_map(|line| {
                    let major = line.split('.').next()?;
                    match major {
                        "18" => Some("Visual Studio 18 2026"),
                        "17" => Some("Visual Studio 17 2022"),
                        "16" => Some("Visual Studio 16 2019"),
                        "15" => Some("Visual Studio 15 2017"),
                        _ => None,
                    }
                })
            })
        })
}

fn msvc_cmake_arch(target_arch: &str) -> Option<&'static str> {
    match target_arch {
        "x86_64" => Some("x64"),
        "x86" => Some("Win32"),
        "aarch64" => Some("ARM64"),
        _ => None,
    }
}

fn run(command: &mut Command, context: &str) {
    let status = command.status().unwrap_or_else(|err| {
        panic!("failed to {context}: {err}", context = context, err = err)
    });
    if !status.success() {
        panic!(
            "{context} failed with status {status}",
            context = context,
            status = status
        );
    }
}

fn built_library_dir(
    build_dir: &Path,
    target_env: &str,
    build_type: &str,
) -> PathBuf {
    let mut candidates = Vec::new();
    if target_env == "msvc" {
        candidates
            .push(build_dir.join("blst").join(build_type).join("blst.lib"));
        candidates.push(build_dir.join("blst").join("blst.lib"));
    } else {
        candidates.push(build_dir.join("blst").join("libblst.a"));
        candidates
            .push(build_dir.join("blst").join(build_type).join("libblst.a"));
    }

    let found: Vec<_> = candidates
        .iter()
        .filter(|candidate| candidate.is_file())
        .collect();
    if found.len() != 1 {
        let expected = candidates
            .iter()
            .map(|path| path.display().to_string())
            .collect::<Vec<_>>()
            .join(", ");
        panic!(
            "expected exactly one built CatBLST archive, found {}; candidates: {}",
            found.len(),
            expected
        );
    }
    found[0]
        .parent()
        .expect("built library must have a parent directory")
        .to_path_buf()
}

fn emit_native_link(lib_dir: &Path, target_env: &str, target_os: &str) {
    println!("cargo:rustc-link-search=native={}", lib_dir.display());
    println!("cargo:rustc-link-lib=static=blst");
    if target_env != "msvc" {
        let cpp_runtime = if target_os == "macos" || target_os == "ios" {
            "c++"
        } else {
            "stdc++"
        };
        println!("cargo:rustc-link-lib={cpp_runtime}");
    }
}

fn emit_fuzz_sanitizer_runtime(target_env: &str, target_os: &str) {
    let Ok(sanitizer) = env::var("BLST_FUZZ_NATIVE_SANITIZER") else {
        return;
    };
    match sanitizer.as_str() {
        "address" => {}
        "undefined" if target_env != "msvc" && target_os == "linux" => {
            println!("cargo:rustc-link-lib=dylib=ubsan");
        }
        "undefined" => panic!(
            "BLST_FUZZ_NATIVE_SANITIZER=undefined is supported only on Linux"
        ),
        _ => panic!(
            "BLST_FUZZ_NATIVE_SANITIZER must be address or undefined"
        ),
    }
}

fn emit_metadata(blst_base_dir: &Path, config_include: &Path) {
    println!(
        "cargo:BINDINGS={}",
        blst_base_dir.join("bindings").to_string_lossy()
    );
    println!(
        "cargo:C_SRC={}",
        blst_base_dir.join("blst").to_string_lossy()
    );
    println!("cargo:CONFIG_INCLUDE={}", config_include.to_string_lossy());
}

fn link_prebuilt(
    library: &Path,
    config_include: &Path,
    target_env: &str,
    target_os: &str,
    blst_base_dir: &Path,
) {
    let expected_name = if target_env == "msvc" {
        "blst.lib"
    } else {
        "libblst.a"
    };
    if library.file_name().and_then(|name| name.to_str()) != Some(expected_name)
    {
        panic!(
            "BLST_PREBUILT_LIB must name {expected_name} for this target: {}",
            library.display()
        );
    }
    if !library.is_file() {
        panic!("BLST_PREBUILT_LIB is not a file: {}", library.display());
    }
    if !config_include.join("config.h").is_file() {
        panic!(
            "BLST_PREBUILT_CONFIG_INCLUDE does not contain config.h: {}",
            config_include.display()
        );
    }
    println!("cargo:rerun-if-changed={}", library.display());
    println!("cargo:rerun-if-changed={}", config_include.display());
    emit_native_link(
        library
            .parent()
            .expect("prebuilt library must have a parent directory"),
        target_env,
        target_os,
    );
    emit_fuzz_sanitizer_runtime(target_env, target_os);
    emit_metadata(blst_base_dir, config_include);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optimized_profiles_select_release_native_code() {
        assert_eq!(native_build_type("0"), "Debug");
        for level in ["1", "2", "3", "s", "z"] {
            assert_eq!(native_build_type(level), "Release");
        }
    }

    #[test]
    fn avx_is_x86_64_only() {
        assert!(validate_build_policy("x86_64", false, true, false).is_ok());
        assert!(validate_build_policy("aarch64", false, true, false).is_err());
    }

    #[test]
    fn conflicting_feature_sets_fail_closed() {
        assert!(validate_build_policy("x86_64", true, true, false).is_err());
        assert!(validate_build_policy("x86_64", true, false, true).is_err());
    }

    #[test]
    fn cross_builds_require_an_existing_toolchain_file() {
        let existing_file = std::env::current_exe().unwrap();
        assert!(validate_cross_toolchain(
            "x86_64-pc-windows-msvc",
            "x86_64-pc-windows-msvc",
            None
        )
        .is_ok());
        assert!(validate_cross_toolchain(
            "x86_64-pc-windows-msvc",
            "aarch64-apple-darwin",
            None
        )
        .is_err());
        assert!(validate_cross_toolchain(
            "x86_64-pc-windows-msvc",
            "aarch64-apple-darwin",
            Some(Path::new("definitely-not-a-toolchain.cmake"))
        )
        .is_err());
        assert!(validate_cross_toolchain(
            "x86_64-pc-windows-msvc",
            "aarch64-apple-darwin",
            Some(existing_file.as_path())
        )
        .is_ok());
    }
}

// Copyright Category Labs. Apache-2.0.
'use strict';
const { spawnSync } = require('node:child_process');
const { createHash } = require('node:crypto');
const [native, wasm] = process.argv.slice(2);
if (!native || !wasm) throw new Error('Usage: node compare-vectors.cjs native-executable wasm-vectors.js');
function run(command, args) {
    const result = spawnSync(command, args, { maxBuffer: 1024 * 1024, timeout: 60000 });
    if (result.error) throw result.error;
    if (result.status !== 0) throw new Error(`${command} failed: ${result.stderr}`);
    return result.stdout.toString('utf8').replaceAll('\r\n', '\n');
}
const expected = run(native, []);
const actual = run(process.execPath, [wasm]);
if (expected !== actual) throw new Error('Native and WebAssembly canonical vectors differ');
console.log(`Native/WASM vectors match: SHA-256 ${createHash('sha256').update(actual).digest('hex')}`);

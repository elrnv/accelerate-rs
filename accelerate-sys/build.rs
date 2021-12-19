extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn build_capi() {
    //clang -c -Wall -Wno-nullability-completeness -o lib.o -I ../src ../src/api.c && ar rc ../libaccelerate_sparse.a lib.o
    cc::Build::new()
        .file("capi/src/api.c")
        .include("capi/src")
        .flag("-Wno-nullability-completeness")
        .compile("libaccelerate_sparse.a")
}

fn main() {
    build_capi();
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    let sdk_root = PathBuf::from("/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk");

    let mut framework_path = sdk_root.clone();
    framework_path.push("System/Library/Frameworks");
    let mut include_path = framework_path.clone();
    include_path.push("Accelerate.framework/Headers");

    println!("cargo:rustc-link-search=native={}", out_path.display());
    println!("cargo:rustc-link-lib=static=accelerate_sparse");

    println!("cargo:rustc-link-lib=framework=Accelerate");

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=capi");

    let bindings = bindgen::Builder::default()
        .header("capi/src/api.h")
        //.detect_include_paths(true)
        .clang_arg(&format!("-F{}", framework_path.display()))
        .clang_arg(&format!("-I{}", include_path.display()))
        .allowlist_function("Sparse.*")
        .allowlist_type("Sparse.*")
        .allowlist_var("Sparse.*")
        .generate_comments(false) // Currently a bug in parsing causes issues.
        .layout_tests(false) // Currently causing UB (https://github.com/rust-lang/rust-bindgen/pull/2055)
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

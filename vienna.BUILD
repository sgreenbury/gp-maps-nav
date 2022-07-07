# https://github.com/bazelbuild/rules_cc/blob/main/examples/my_c_archive/my_c_archive.bzl
# https://docs.bazel.build/versions/main/cpp-use-cases.html
cc_library(
    name = "RNA",
    srcs = glob(["lib/*.a"]),
    hdrs = glob(["include/ViennaRNA/*.h"]),
    copts = ["-Iinclude/ViennaRNA/", "-Llib/"],
    # , "-lm", "-lRNA"],
    linkstatic = 1,
    visibility = ["//visibility:public"],
)

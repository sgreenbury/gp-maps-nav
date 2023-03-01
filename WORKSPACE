# Import from symlink of installed static ViennaRNA library
new_local_repository(
    name = "vienna_rna",
    path = "ViennaRNA/",
    build_file = "vienna.BUILD",
)

# Import http_archive
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# absl
http_archive(
  name = "com_google_absl",
  urls = ["https://github.com/abseil/abseil-cpp/archive/215105818dfde3174fe799600bb0f3cae233d0bf.zip"],
  strip_prefix = "abseil-cpp-215105818dfde3174fe799600bb0f3cae233d0bf",
  sha256 = "b4e20d9e752a75c10636675691b1e9c2698e0764cb404987d0ffa77223041c19"
)

# C++ rules for Bazel
http_archive(
  name = "rules_cc",
  urls = ["https://github.com/bazelbuild/rules_cc/archive/262ebec3c2296296526740db4aefce68c80de7fa.zip"],
  strip_prefix = "rules_cc-262ebec3c2296296526740db4aefce68c80de7fa",
  sha256 = "9a446e9dd9c1bb180c86977a8dc1e9e659550ae732ae58bd2e8fd51e15b2c91d",
)

# GoogleTest
http_archive(
  name = "com_google_googletest",
  urls = ["https://github.com/google/googletest/archive/58d77fa8070e8cec2dc1ed015d66b454c8d78850.zip"],
  strip_prefix = "googletest-58d77fa8070e8cec2dc1ed015d66b454c8d78850",
  sha256 = "ab78fa3f912d44d38b785ec011a25f26512aaedc5291f51f3807c592b506d33a"
)

# Benchmark
http_archive(
    name = "com_github_google_benchmark",
    urls = ["https://github.com/google/benchmark/archive/bf585a2789e30585b4e3ce6baf11ef2750b54677.zip"],
    strip_prefix = "benchmark-bf585a2789e30585b4e3ce6baf11ef2750b54677",
    sha256 = "2a778d821997df7d8646c9c59b8edb9a573a6e04c534c01892a40aa524a7b68c",
)

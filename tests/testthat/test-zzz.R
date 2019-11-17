test_that("`fetch_java_version()` works", {
  expect_is(fetch_java_version(), "character")
})

test_that("`check_java_version() works", {
  expect_null(check_java_version())
})

test_that("`check_java_version() works with 1.8", {
  expect_null(check_java_version(
    c("java version \"1.8.0_144\"",
      "Java(TM) SE Runtime Environment (build 1.8.0_000-000)",
      "Java HotSpot(TM) 64-Bit Server VM (build 00.000-000, mixed mode)")))
})

test_that("`check_java_version() fails with 1.7", {
  expect_error(check_java_version(
    c("java version \"1.7.0\"",
      "Java(TM) SE Runtime Environment (build 1.7.0_000-000)",
      "Java HotSpot(TM) 64-Bit Server VM (build 00.000-000, mixed mode)")))
})

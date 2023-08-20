## Tests for functions related to java version check - Aug 2023

test_that("`fetch_java_version()` works as expected", {
    version_vec <- c("java version \"13.0.1\" 2019-10-15", "Java(TM) SE Runtime Environment (build 13.0.1+9)",
        "Java HotSpot(TM) 64-Bit Server VM (build 13.0.1+9, mixed mode, sharing)")

    mockery::stub(fetch_java_version, "Sys.getenv", "/path/to/java/home")
    mockery::stub(fetch_java_version, "file.exists", TRUE)
    mockery::stub(fetch_java_version, "system2", version_vec)

    # unix
    mockery::stub(fetch_java_version, "identical", FALSE)
    expect_equal(fetch_java_version(), version_vec)
    # windows
    mockery::stub(fetch_java_version, "identical", TRUE)
    expect_equal(fetch_java_version(), version_vec)

    mockery::stub(fetch_java_version, "system2", c())
    expect_error(fetch_java_version())

    mockery::stub(fetch_java_version, "file.exists", FALSE)
    expect_error(fetch_java_version())

    mockery::stub(fetch_java_version, "Sys.getenv", NA)
    mockery::stub(fetch_java_version, "Sys.which", "path/to/java")
    mockery::stub(fetch_java_version, "system2", version_vec)
    expect_equal(fetch_java_version(), version_vec)
})

test_that("`check_java_version()` works", {
    expect_null(check_java_version())
})

test_that("`check_java_version()` raises parsing error", {
    expect_error(check_java_version(c("version 1.8", "version 1.7")), "Java version detected but couldn't parse version from ")
    expect_error(check_java_version("version XXXX"), "Java version detected but couldn't parse version from: ")
})

test_that("`check_java_version()` works with 1.8", {
    expect_null(check_java_version(c("java version \"1.8.0_144\"", "Java(TM) SE Runtime Environment (build 1.8.0_000-000)",
        "Java HotSpot(TM) 64-Bit Server VM (build 00.000-000, mixed mode)")))
})

test_that("`check_java_version()` works with 14", {
    expect_null(check_java_version(c("java version \"14\" 2020-03-17", "Java(TM) SE Runtime Environment (build 14+36-1461)",
        "Java HotSpot(TM) 64-Bit Server VM (build 14+36-1461, mixed mode, sharing)")))
})

test_that("`check_java_version()` fails with 1.7", {
    expect_error(check_java_version(c("java version \"1.7.0\"", "Java(TM) SE Runtime Environment (build 1.7.0_000-000)",
        "Java HotSpot(TM) 64-Bit Server VM (build 00.000-000, mixed mode)")))
})

find_package(OpenSSL REQUIRED)
add_definitions(
  -DCURVE_BN128
)

add_subdirectory(src)

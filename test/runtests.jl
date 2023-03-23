using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

const is_TRAVIS = haskey(ENV, "TRAVIS")

const is_CI = haskey(ENV, "CI")

# Currently verified by Downstream tests
@test true

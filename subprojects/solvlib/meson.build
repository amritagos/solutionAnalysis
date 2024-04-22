project('solvlib', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20'])

_deps = [] # Dependencies

# Dependencies
# For clang
_deps += [declare_dependency(link_args: '-lstdc++')]

# fmt
fmt_dep = dependency('fmt', required: true)
_deps += fmt_dep

inc = include_directories('include')
solvlib = library('simple',
  'add.cpp',
  include_directories : inc,
  dependencies: _deps,
  install : true)

solvlib_dep = declare_dependency(include_directories : inc,
  link_with : solvlib, dependencies: _deps)

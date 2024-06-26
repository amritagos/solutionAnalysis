project('soluanalysis', 'cpp', 
  default_options : ['warning_level=3', 'cpp_std=c++20'])

# Add C++ compiler options
cpp_args = ['-Wall', '-Wextra', '-Wpedantic']
_args = [cpp_args] # Extra arguments
_deps = [] # Dependencies

# ---------------------- Library Dependencies

solvlib_dep = dependency('solvlib', fallback : ['solvlib', 'solvlib_dep'])
_deps += solvlib_dep

graphlib_dep = dependency('graph_lib', fallback : ['graph_lib', 'graphlib_dep'])
_deps += graphlib_dep

# ---------------------- Bindings

py_mod = import('python')
py = py_mod.find_installation(pure: false)
pyb11_dep = [
  py.dependency(),
  dependency('pybind11')
]
_deps += [pyb11_dep]

# for the bindings
py.extension_module(
  'solvlib',
  sources : [
    'python_bindings/bindings.cpp'
  ],
  dependencies: _deps,
  cpp_args : _args,
  install: true,
  subdir: 'soluanalysis/'
)

py.extension_module(
  'graphlib',
  sources : [
    'python_bindings/bindings.cpp'
  ],
  dependencies: _deps,
  cpp_args : _args,
  install: true,
  subdir: 'soluanalysis/'
)

# soluanalysis, main package
py.install_sources([
    'soluanalysis/__init__.py',
  ],
  pure: false, # install next to compiled extension
  subdir: 'soluanalysis'
)
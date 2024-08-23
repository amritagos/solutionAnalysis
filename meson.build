project('soluanalysis', 'cpp', 
  default_options : ['warning_level=3', 'cpp_std=c++20'])

# Add C++ compiler options
cpp_args = ['-Wall', '-Wextra', '-Wpedantic']
_args = [cpp_args] # Extra arguments
_deps = [] # Dependencies

# ---------------------- Library Dependencies

# james has graph_lib as a dependency already
james_dep = dependency('james', fallback : ['james', 'james_dep'])
_deps += james_dep

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
  'james',
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

# IO and extras
py.install_sources([
    'soluanalysis/io.py',
    'soluanalysis/ion_pairs.py',
    'soluanalysis/misc.py'
  ],
  pure: false,
  subdir: 'soluanalysis'
)
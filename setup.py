from distutils.core import setup, Extension

pygem_module=Extension('pygem',
                       library_dirs=['/usr/local/lib'],
                       libraries=['gem'],
                       include_dirs=['/usr/local/include'],
                       sources=['pygem.c'])
setup(name='pygem',version='1.0',description='Python extension of the libgem library',ext_modules=[pygem_module])

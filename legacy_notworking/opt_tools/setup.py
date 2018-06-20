from distutils.core import setup, Extension


module1 = Extension('libOPTcalcs',
                    define_macros = [('MAJOR_VERSION', '1'),('MINOR_VERSION', '0')],
                    sources = ['libOPT_calcs.c'])

setup (name = 'libOPTcalcs',
              version = '1.0',
              description = 'This is the calculations for libOPT',
              author = 'Johnathon R. Walls',
              author_email = 'johnathon.walls@gmail.com',
              long_description = '''
              This is just the calculations for libOPT.
              ''',
              ext_modules = [module1])

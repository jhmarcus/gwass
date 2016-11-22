from setuptools import setup

setup(name='gwass',
      version='0.1',
      description='',
      url='https://github.com/jhmarcus/gwass',
      author='Joseph Marcus',
      author_email='jhmarcus@uchicago.edu',
      license='MIT',
      packages=['gwass'],
      install_requires=[
          'numpy',
          'pandas',
          'cython',
          'pysam'
      ],
      zip_safe=False)

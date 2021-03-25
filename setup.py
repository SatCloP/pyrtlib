from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='pyrtlib',
    version='1.0.0',
    packages=['pyrtlib'],
    install_requires=required,
    url='',
    license='MIT License',
    author='slarosa',
    author_email='salvatore.larosa@imaa.cnr.it',
    description='pyrtlib - Radiative Transfer library'
)

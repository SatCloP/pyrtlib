from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='pyrtlib',
    version='1.0.0',
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.6',
    install_requires=required,
    url='',
    license='MIT License',
    author='slarosa',
    author_email='salvatore.larosa@imaa.cnr.it',
    description='pyrtlib - Radiative Transfer library'
)

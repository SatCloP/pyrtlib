from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='pyrtlib',
    version='1.3.0',
    include_package_data=True,
    package_dir={'': 'pyrtlib'},
    python_requires='>=3.6',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=required,
    url='',
    license='MIT License',
    author='slarosa',
    author_email='salvatore.larosa@imaa.cnr.it',
    description='pyrtlib - Radiative Transfer library',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering'
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    packages=find_packages('pyrtlib')
)

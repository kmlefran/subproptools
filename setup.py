import setuptools

setuptools.setup(
    name='subproptools',
    version='0.0.3',
    description='Extract group properties from .sum file and rotate substituent geometries',
    long_description=open('README.md').read(),
    url='https://github.com/kmlefran/subproptools',
    author='Kevin M. Lefrancois-Gagnon',
    install_requires=['numpy','pandas'],
    author_email='kgagnon@lakeheadu.ca',
    packages=setuptools.find_packages(),
    classifiers = [
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.11',
    ],
    zip_safe=False,
    include_package_data=True,
    package_data={'': ['test_data/*.sum']},
)
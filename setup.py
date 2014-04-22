from distutils.core import setup

setup(
        name='fasta2tab',
        version='0.1.0',
        author='J. D. A. Gilliland',
        author_email='jdagilliland@gmail.com',
        py_modules=['fasta2tab'],
        scripts=['bin/fasta2tab'],
        install_requires=[
            'BioPython >= 1.63',
            ],
        license='LICENSE.txt',
        description='Utilities to help with visualizing PHYLIP trees',
        long_description=open('README.md').read(),
        )



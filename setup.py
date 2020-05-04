import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyDSM-afeborgeaud",
    version='0.0.1',
    author='Anselme Borgeaud',
    author_email='aborgeaud@gmail.com',
    description='Python wrapper for DSM',
    long_description = long_description,
    long_description_content_type='text/markdown',
    url='',
    packages=['pyDSM', ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7'
)
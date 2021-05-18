from setuptools import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
    long_description_content_type = "text/x-rst"
except(IOError, ImportError):
    long_description = open('README.md').read()
    long_description_content_type = 'text/markdown'

setup(
    name='zCurve',
    version='0.0.3',
    description='zCurve maps multidimensional data to one dimension while preserving locality of the data points.',
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    url='https://github.com/rmrschub/pyMorton',
    author='René Schubotz',
    author_email='rene.schubotz@dfki.de',
    license='CC BY-NC-SA 4.0',
    packages=['zCurve'],
    install_requires=['typing',
                      'gmpy2'],

    classifiers=[ ],
)
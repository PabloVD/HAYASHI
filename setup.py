# Following https://www.kdnuggets.com/pip-install-you-a-beginners-guide-to-creating-your-python-library

from setuptools import find_packages, setup

with open('README.md') as f:
    long_description = f.read()

setup(
    name='hayashi',
    packages=find_packages(),
    version='1.1.3',
    description='Python library for computing the number of absorption features of the 21 cm forest in a semianalytic formalism.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Pablo Villanueva-Domingo',
    author_email='pablo.villanueva.domingo@gmail.com',
    install_requires=["numpy",
                      "scipy",
                      "tqdm",
                      "colossus"],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # License type
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    url='https://github.com/PabloVD/HAYASHI',
)

from setuptools import setup, find_packages

setup(
    name='code',
    version='1.0.0',    
    description='The pymCADRE tool is an advanced re-implementation of the metabolic Context-specificity Assessed by Deterministic Reaction Evaluation (mCADRE) algorithm in Python. It constructs tissue-specific metabolic models by leveraging gene expression data and literature-based evidence, along with network topology information.',
    url='https://github.com/draeger-lab/pymCADRE',
    author='Nantia Leonidou',
    author_email='nantia.leonidou@uni-tuebingen.de',
    license=' GPL-3.0',
    packages=find_packages(include=['rank', 'check','prune']),
    install_requires=['pandas',
                      'numpy',
                      'cobrapy',
                      'requests'
                      'os'
                      ],

    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Topic :: Software Development :: Version Control :: Git",
        "Framework :: Jupyter",
        "Operating System :: MacOS",
        "Operating System :: Unix"
   
    ],
)
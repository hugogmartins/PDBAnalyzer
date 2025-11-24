from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pdb_analyzer",
    version="0.1.0",
    author="Hugo Guimarães Martins, Adenilson Arcanjo, Milena Pirovani e Valentina Lourenço",
    author_email="hugoguimaraesmartins@gmail.com",
    description="Protein structure analysis toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hugogmartins/PDBAnalyzer",
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific :: Bio-Informatics",
        "Topic :: Scientific :: Bio-Chemistry",
    ],
    
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    
    install_requires=[
        "biopython>=1.81",
        "numpy>=1.21.0",
        "matplotlib>=3.5.0",
    ],
    
    entry_points={
        "console_scripts": [
            "protein-analyzer=protein_analyzer.cli:main",
        ],
    },
    
    include_package_data=True,
    package_data={
        "protein_analyzer": ["config/*.yaml", "templates/*.html"],
    },
    
    keywords="bioinformatics protein structure analysis PDB biopython",
)
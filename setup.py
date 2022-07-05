import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FunFun",
    version="0.1.15",
    author="D.V. Krivonos, D.N. Konanov",
    author_email="danil01060106@gmail.com",
    description="ITS-based functional annotator",
    long_description="FunFun",
    long_description_content_type="",
    url="https://github.com/DanilKrivonos/FunFun",
    project_urls={
        "Bug Tracker": "https://github.com/DanilKrivonos/FunFun",
    },
    package_data={
        'mypackage': ['FunFun/data/CONCATENATE_base.tsv'],
	'mypackage': ['FunFun/data/functionality.tsv'],
	'mypackage': ['FunFun/data/ITS1_base.tsv'],
	'mypackage': ['FunFun/data/ITS2_base.tsv'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    include_package_data=True,
    packages=['FunFun', 'FunFun.src', 'FunFun.example', 'FunFun.data'],
    install_requires=[
        'numpy',
        'biopython',
        'pandas',
        'plotly',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'funfun=FunFun.FunFun:main'
        ]
    }
)

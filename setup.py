import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FunFun",
    version="0.1.0",
    author="D.V. Krivonos, D.N. Konanov",
    author_email="danil01060106@gmail.com",
    description="ITS-based functional annotator",
    long_description="FunFun",
    long_description_content_type="",
    url="https://github.com/DanilKrivonos/FunFun",
    project_urls={
        "Bug Tracker": "https://github.com/DanilKrivonos/FunFun",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    include_package_data=True,
    packages=['FunFun.py', 'FunFun.src', 'FunFun.external', 'FunFun.data'],
    install_requires=[
        'numpy',
        'biopython',
        'pandas',
        'plotly'
    ],
    entry_points={
        'console_scripts': [
            'funfun=FunFun.FunFun:main'
        ]
    }
)
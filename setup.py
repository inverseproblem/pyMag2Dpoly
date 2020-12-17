

from setuptools import setup, find_packages

setup(
    version="0.2.0"
    name="Mag2DPoly",
    author="Andrea Zunino, Alessandro Ghirotto",
    description="Forward magnetic anomaly calculation due to two-dimensional polygonal bodies with uniform arbitrary polarization",
    long_description=open('README.md').read(),
    license="LICENSE",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=["numpy"], # "matplotlib"],
    extras_require={
        "dev": [
            #Examples
            "matplotlib",
            # Documentation
            "sphinx",
            "sphinx_rtd_theme",
        ]
    },
    zip_safe=False,
)

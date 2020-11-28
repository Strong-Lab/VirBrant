<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/Strong-Lab/VirBrant">
    <img src="images/logo.jpg" alt="Logo" width="160" height="70">
  </a>

  <h3 align="center">VirBrant</h3>

  <p align="center">
    A viral identification tool using machine learning with nucleotide and protein features
    <br />
    <a href="https://github.com/Strong-Lab/VirBrant/issues">Report Bug</a>
    ·
    <a href="https://github.com/Strong-Lab/VirBrant/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-VirBrant">About VirBrant</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About VirBrant

This project was created to identify viral contigs in metagenomics. Add details about project and images

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

VirBrant requires Python 3 and the following libraries (if installling through pip, libraries are automatically install)
* pandas
* scikit-learn
* biopython


### Installation

VirBrant can only be forked from this repository. Once forked, enter into the files and compile the kmer counting program. 

```
## git download here

cd VirBrant/kmer_counter
make
```

<!-- USAGE EXAMPLES -->
## Usage
VirBrant works as a command line script. Once install via pip, VirBrant the command can be accessed. To get the help screen type:
```
VirBrant -h
```

The paramters of VirBrant are:
* -f: Kraken output file \[required]
* -c: Seqeuncing file to parse \[optional]
* -r: Remove viral elements flag
* -o: Rename output files \[optional]


### Running VirBrant 

#### VirBrant without a sequencing file and renaming the output
```
VirBrant -f Kraken_Output.txt -o Viral_Sequences
```


<!-- ROADMAP -->
## Roadmap

<p align="center">
    Current Version: 0.0.1
</p>

Improvements to be made:
- Reduce feature space to allow for smaller file processing
- Grid search hyper parameters for models
- Build into python package


See the [open issues](https://github.com/othneildrew/Best-README-Template/issues) for a list of proposed features (and known issues).


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Cody Glickman - [@glickman_Cody](https://twitter.com/glickman_cody) - glickman.cody@gmail.com

Project Link: [https://github.com/Strong-Lab/VirBrant](https://github.com/Strong-Lab/VirBrant)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* James Costello
* Michael Strong
* Jo Hendrix





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/Strong-Lab/VirBrant.svg?style=for-the-badge
[contributors-url]: https://github.com/ontributors/Strong-Lab/VirBrant/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/Strong-Lab/VirBrant.svg?style=for-the-badge
[forks-url]: https://github.com/Strong-Lab/VirBrant/network/members
[stars-shield]: https://img.shields.io/github/stars/Strong-Lab/VirBrant.svg?style=for-the-badge
[stars-url]: https://github.com/Strong-Lab/VirBrant/stargazers
[issues-shield]: https://img.shields.io/github/issues/Strong-Lab/VirBrant.svg?style=for-the-badge
[issues-url]: https://github.com/Strong-Lab/VirBrant/issues
[license-shield]: https://img.shields.io/github/license/Strong-Lab/VirBrant.svg?style=for-the-badge
[license-url]: https://github.com/Strong-Lab/VirBrant/LICENSE.txt


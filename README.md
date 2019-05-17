# pep2struc
A tool that help converting peptide-sequences into their 2D structures.

## Installation
This script uses Python 3.7.x. If you don't have Python, I would recommend downloading it from [Anaconda](https://www.continuum.io/downloads).

Copy or clone this package from Github.

Open the Terminal/Command Line and navigate to where you copied the package:

    $ cd path/to/copied/directory

### Linux and MacOS

Install the dependencies by entering:

    $ pip install -r requirements.txt

## Usage

To run from the command-line, just do:

    $ python pep2struc.py

Example: Creating linear peptides

    $ python pep2struc.py -in path_to/peptideSequences.csv -o path_to_output/ -t linear

To list all the parameters needed from the command-line, just do:

    $ python pep2struc.py --help

## Questions and Comments

Feel free to direct any questions or comments to the Issues page of the repository.

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).

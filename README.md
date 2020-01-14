# informative_sites

Python code to gather informative sites from sequence alignment .fasta format file, and write it out to .fasta file.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Python version > 3.4
smallBixTools
argparse
```

### Installing

Clone the repo.
Run it from command line like so:

```
python3 informative_sites.py -in /path/to/input/fasta.fasta -out_dir /path/to/output/ --caseSensitive
```

Example output when used on the example_input.fasta contained with this code:

```
Input alignment is 73 sites long.
After removing conserved sites, the alignment has 5 sites.
Completed.                                                                                                                                                          
Now exiting
```


## Running the tests

Currently no tests exist.


## Built With

* [PyCharm](https://www.jetbrains.com/pycharm/) - The IDE used


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We want to use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/HIVDiversity/informative_sites/tags).

## Authors

* **David Matten** - *Initial work* - [HIV Diversity Group](https://github.com/HIVDiversity/)

See also the list of [contributors](https://github.com/HIVDiversity/informative_sites/graphs/contributors) who participated in this project.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Hat tip to Melissa-Rose for the project inspiration (https://github.com/orgs/HIVDiversity/people/Melissa-Rose)



# rustracking

`rustracking` is a personal project developed with the intent of learning Rust. I started working on it after having read [the Book](https://doc.rust-lang.org/stable/book/) and wanting to create a project in a (fun) topic I am familiar with.
It's final goal is to perform a simplified version of track reconstruction in particle physics experiments.
Everything should be implemented from scratch as much as possible, using only the standard library.

## Progress
Some of the things which will have to be implemented

- [x] Matrix
- [x] Module + IO
- [x] Hits + IO
- [x] \*Event generation
- [x] Clustering
- [ ] Seeding
- [ ] Visual representation

\* done but not yet satisfied with result.

## Contributing

If you happen to stumble upon this project, feel free help out, keeping in mind that it is a learning exercise.

## Progress pictures

The box detector with randomly generated rays, resulting truth hits (green), activated cells, and clustered 3d points (black):

<img src="https://github.com/guilhermeAlmeida1/rustracking/blob/d91483292871ee5a8ecc09f806433a4e734ea324/data/randomHits/3d.svg" width="560" />

A _very cool_ block of code, where I decided to process the output of the CCL algorithm in a functional manner:
[clustering.rs#L85-L106](https://github.com/guilhermeAlmeida1/rustracking/blob/d91483292871ee5a8ecc09f806433a4e734ea324/src/clustering.rs#L85-L106)

Magnetic field integration into event generator:

<img src="https://github.com/guilhermeAlmeida1/rustracking/blob/87b12a0fb1d7a16f6c6ebe02da48e303a0a6dfe7/data/magfield/3d.svg" width="560" />

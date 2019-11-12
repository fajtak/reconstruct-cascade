# An example of a BARS standalone app

This repository hosts an example of a BARS app - a compiled program that we advise should be used instead of ROOT macros.

## Requirements

This app only requires a BARS installation. 
The exact steps are described [here](https://bars-docs.jinr.ru/install).

## Installation and compilation 

1. Clone this repository to a directory named as this app

```
git clone git@git.jinr.ru/Baikal/bars-app-standalone-example.git my-bars-app
```

2. Create a build directory.

```
cd my-bars-app
mkdir build
```

3. Create a Makefile and compile

```
cd build
cmake ..
make
```

4. Install it

```
make install
```

This will put an output binary into `~/.local/bin/bars`. Make sure `~/.local/bin` is  in your `$PATH`!.


4. Start your app!

```
my-bars-app
```

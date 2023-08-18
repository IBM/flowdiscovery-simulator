# Developer guidelines

These guidelines are a sure way of helping you develop the code while maintaining consistency with
the current style, structure and formatting choices.

***

## Git branching model

### [TL;DR](https://en.wiktionary.org/wiki/TLDR) version

This project uses a [*shared repository model*](https://help.github.com/articles/about-collaborative-development-models/),
in which feature branches are created and merged into `master` via [Pull Requests](https://help.github.com/articles/about-pull-requests/) (PR).

1. Pick an issue from the [list](https://github.com/IBM/flowdiscovery-simulator/issues).
1. Create a local `feature_branch` named after the issue.
1. Code, commit & push to `remotes/origin/feature_branch`.
1. Build and run the regression tests (`./bin/flow-simulator-test.x`).
1. If the build succeeds, issue a Pull Request towards `master`.
1. After the PR is accepted, remove the local `feature_branch`.
1. Update your local `master` and restart from the beginning.

### And remember to, at least once a day

* Pull the latest content of `master` to you local repository and `git merge master` from your `feature_branch`.
* Run regression tests (`./bin/flow-simulator-test.x`) locally.
* Run linter checks.

### Commit messages

Please give meaningful titles to your commit messages and write `Issue #X` in the second line so
that the commit is linked to the issue it addresses.
If necessary add a detailed description of the line of thought behind the proposed implementation
for bookkeeping purposes.

**Example:**

```Text
Added method to calculate The Answer to Life, the Universe and Everything
Issue #42

There is a theory which states that if ever anyone discovers exactly what the Universe is for and
why it is here, it will instantly disappear and be replaced by something even more bizarre and
inexplicable.

There is another theory mentioned, which states that this has already happened.
```

### Pull request messages

Please give meaningful titles to your Pull Requests and write `Closes #X` in the body so that the
issue is automatically closed when the PR is accepted.
If necessary add a detailed description of the line of thought behind the proposed implementation
for bookkeeping purposes.

### [TMI](https://en.wiktionary.org/wiki/too_much_information#English) version

The `master` branch is protected and can only be altered via Pull Requests from feature branches.
The `master` branch should **always** be deployable at any point in time.

Typically, one feature branch will be created for each issue.
However, for chained issues, it is preferable to address them all in the same branch, followed by a
single PR containing the entire change set.
Feature branches should follow the name convention represented by `meaningful_keywords-X`.
The choice for `meaningful_keywords` should refer to the issue it is supposed to solve.
Likewise, `X` is the number of the target issue.
When dealing with chained issues (`X`, `Y` and `Z`), `meaningful_keywords-X-Y-Z` is also possible.

## Coding style

The only rule is **to be consistent to the style around you**: headers, documentation, include
order, whitespace, alignment, comments, naming of files, variables and functions.
No one should be able to guess who wrote the file based on how it looks like.

### C++ coding style

All **C++** source files in this project follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
This particular style was chosen simply because it offers a tool ([cpplint](https://github.com/cpplint/cpplint)) to check for
compliance.

```Shell
cpplint --recursive --quiet --linelength=100 --extensions=h,cc --filter=-runtime/references,-whitespace/empty_loop_body src/ test/
```

### Python coding style

All **Python** source files in the project follow the [Flake8 Style Guide](http://flake8.pycqa.org/en/latest/manpage.html).

```Shell
flake8 --max-line-length=100 util/ api/
```

## Building with `cmake`

Please configure `cmake` to perform **out-of-source** builds, into the `build` directory.
Besides the standard build targets, `cmake` will generate the following custom targets:

* `cpplint`: Runs `cpplint` on the **C++** files.
* `flake8`: Runs `flake8` on the **Python** files.
* `docs`: Runs `doxygen` and builds documentation at `docs/html/index.html`.
* `clean-docs`: Cleans the content of `docs/`.
* `clean-results`: Cleans the content of `results/`.
* `clean-cmake`: Cleans the content of `build/`.
* `clean-all`: Performs `clean`, `clean-docs`, `clean-results` and `clean-cmake`.

## Documentation

Please document files, classes, functions, methods and members using [Doxygen](https://www.doxygen.nl/index.html).
**Be consistent to the documentation around you**.

* **.cc and .h files** should contain a header with the `\file`, `\brief`, `\authors`, `\copyright`
and `\date` entries.
The long description (or `\details`) should list the classes (`*.h`) or methods (`*.cc`) therein.
Header guards (`*.h`) should represent the file name.

* **Classes** should contain a header with `\class`, `\brief` and a one-line description of what
the class does.

* **Method declarations** should contain a one-line description of what the method does.

* **Method implementations** should contain a detailed description of what the method does and how
it does it, specially, when there are mathematical equations involved.
External references (articles and/or Wikipedia) are welcome.
The input (`\param[in]`) and output (`\param[out]`) parameters as well as the return values
(`\return`) are also required.

* **Member variables** should contain a `\brief` description of what the variable means and a
detailed description of where it is used and for what.

## Continuous integration

After each `git push`, Travis CI will

* lint the source code,
* build the Kubernetes container,
* run regression tests,

**Tip**: If you do not want to run a build for a particular commit for any reason, add `[ci skip]`
or `[skip ci]` to the git commit message.

## Choosing Open Source libraries

This project makes use of open source libraries, as described in the [README.md](README.md) file.

## Source code editor

The project contains an `.editorconfig` file to enforce the indent type (**2 spaces**) and line
length (**100 columns**).
Please pick an editor that supports [EditorConfig](http://editorconfig.org/#download) either
natively or via an extension.

### Visual Studio Code

**Pros:**

* Built-in debugger.
* Faster when starting and building projects.
* Better C++ integration with [cpptools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools).
* Integrated `cpplint` support with [linter-cpplint](https://marketplace.visualstudio.com/items?itemName=mine.cpplint).
* Integrated `flake8` support with [python](https://marketplace.visualstudio.com/items?itemName=donjayamanne.python).
* Good control over `cmake` features with [cmake-tools](https://marketplace.visualstudio.com/items?itemName=vector-of-bool.cmake-tools).
* Support for EditorConfig with [EditorConfig](https://marketplace.visualstudio.com/items?itemName=EditorConfig.EditorConfig).

**Cons:**

* Harder to configure extensions.

**Code --> Preferences --> Settings:**

```JSON
// settings.json
"editor.rulers": [100],
"cpplint.cpplintPath": "cpplint",       # which cpplint
"cpplint.lintMode": "workspace",
"cpplint.lineLength": 100,
"cpplint.verbose": 1,
"cpplint.filters": ["-runtime/references", "-whitespace/empty_loop_body"],
"python.testing.pytestEnabled": false,
"python.testing.unittestEnabled": true,
"python.testing.nosetestsEnabled": false,
"python.linting.pylintEnabled": false,
"python.linting.flake8Enabled": true,
"python.linting.flake8Args": ["--max-line-length=100"],
"C_Cpp.default.cppStandard": "c++14",
"cmake.generator": "Ninja",
"cmake.configureArgs": ["-DCMAKE_C_COMPILER=gcc-11", "-DCMAKE_CXX_COMPILER=g++-11"],
"git.confirmSync": false
```

## Development environment

In case you are not using the same OS as the production environment, i.e. the latest **Fedora**, a
virtual machine can be created using the `Vagrantfile` contained in the root project directory.

```Shell
# On the host machine:
(host) : flow-simulator > vagrant up
(host) : flow-simulator > vagrant ssh

# On the virtual machine:
(VM) : ~/flow-simulator > cmake-generate  # generates Ninja build files
(VM) : ~/flow-simulator > cmake-build     # compiles flow-simulator.x
(VM) : ~/flow-simulator > cmake-build-all # compiles all binaries
(VM) : ~/flow-simulator > cmake-clean     # cleans object files
(VM) : ~/flow-simulator > cmake-clean-all # cleans object and CMake files
```

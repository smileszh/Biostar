# (PART) 模块化脚本 {-}



# makefile

## 简介
官方教程：

https://www.gnu.org/software/make/manual/make.html

在生信上游分析中，我们不可能每次都手动输入我们的命令，通常我们会将我们的命令写入到一个脚本文件中，然后通过运行脚本文件来执行我们的命令。但是，当我们工作有些变化时，我们就需要做出一些改变，包括增加或删除，这样就会增加我们的工作量和代码复用性。因此，我们可以使用 makefile 来模块化我们的分析流程。

Makefile 是一种用于自动化构建过程的文件，通常用于 Unix 和类 Unix 系统上。它由一个名为 `make` 的工具读取和执行。Makefile 定义了一系列的规则和依赖关系，用于指导如何编译和链接程序。

在生信分析中使用 makefile 可以帮助我们自动化分析流程，我们不需要特别高级的学习。
下面是一个简单的 makefile 示例：

首先我们创建一个文件，命名为`Makefile`，然后在文件中写入一些内容。

``` makefile
.RECIPEPREFIX = >

# 命令 1 
foo:
> echo Hello world!

# 命令 2
bar:
> echo Hello world!
> echo Hello Everyone!
```
在这个 makefile 中，我们定义了两个目标 `foo` 和 `bar`，分别对应两个命令。我们可以通过 `make foo` 和 `make bar` 来执行这两个命令。


``` bash
$ make -f Makefile foo
echo Hello world!
Hello world!
```


``` bash
$ make -f Makefile bar
echo Hello world!
Hello world!
echo Hello Everyone!
Hello Everyone!
```
当你的文件名是 `Makefile` 时，你可以直接使用 `make` 命令，而不需要 `-f` 参数。

通常 makefile 会将命令和结果同时输出到终端，如果我们只想输出结果，可以将命令前加一个@符号。比如这样：

``` makefile
.RECIPEPREFIX = >

# 命令 1 
foo:
> @echo Hello world!
```

和 bash 脚本一样，makefile 也可以使用变量，比如：

``` makefile
.RECIPEPREFIX = >

NAME = xiaoming

hello:
> @echo Hello ${NAME}
```
输出结果为：

``` bash
$ make -f anyfile hello
Hello xiaoming
```

`make` 的开发目的是跟踪目标之间的相互依赖关系，并在任何文件发生更改时，仅执行必要的步骤。例如，如果你已经有一个基因组索引，`make` 命令将跳过索引步骤。

``` makefile
# Set the prefix from tabs to >
.RECIPEPREFIX = >

counts.txt:
> echo 100 > counts.txt

names.txt:
> echo Joe > names.txt

results.txt: counts.txt names.txt
> cat counts.txt names.txt > results.txt
```
在上面的 Makefile 中，输入 `make results.txt` 将首先生成 `counts.txt` 和 `names.txt`，如果它们不存在的话。然后，它将从这两个文件生成 `results.txt`。如果你再次运行 `make`，它将不会执行任何操作，因为所有文件已经存在。

## makefile 的小技巧

### 默认设置
在 Makefile 中添加上这些设置后可以让 makefile 更加强大：

``` makefile
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory
```
`.DELETE_ON_ERROR:`

这个特殊目标告诉 `make` 在命令执行失败时删除生成的目标文件，以防止生成不完整或损坏的文件。
`SHELL := bash`

这行代码指定 `make` 使用 `bash` 作为默认的 shell 来执行命令。默认情况下，`make` 使用 `/bin/sh`，但可以通过这种方式指定使用其他 shell。
`.ONESHELL:`

这个特殊目标指示 `make` 将单个目标的所有命令在同一个 shell 会话中执行，而不是为每个命令启动一个新的 shell。这对于需要在多个命令之间共享环境或状态的情况很有用。
`.SHELLFLAGS := -eu -o pipefail -c`

这些标志用于配置 `bash` 的行为：
`-e`：如果任何命令失败（返回非零状态），则 `bash` 退出。
`-u`：使用未定义的变量时，`bash` 会报错并退出。
`-o pipefail`：如果管道中的任何命令失败，整个管道返回失败状态。
`-c`：从字符串中读取命令并执行。
`MAKEFLAGS += --warn-undefined-variables --no-print-directory`

`--warn-undefined-variables`：警告使用未定义的变量，以帮助捕获拼写错误或遗漏的变量定义。
`--no-print-directory`：禁用 `make` 在递归调用时打印目录信息，这可以使输出更简洁。


### “试运行”模式
在 bash 命令中，命令对不对，我们可以运行一下看看，如果不对，我们可以再次运行。但是在 makefile 中，我们不希望这样，我们希望一次就对，所以我们可以使用 `-n` 参数来进行“试运行”模式。这样 makefile 会输出将要执行的命令，但是不会真正执行。这样可以帮助我们检查命令是否正确。


``` bash
$ make hello
Hello xiaoming

$ make hello -n
echo Hello xiaoming
```

### 变量命名
在 Makefile 中可以通过 `?=` 这样的方式预设变量：

``` makefile
.RECIPEPREFIX = >

NAME ?= xiaoming

hello:
> @echo Hello ${NAME}
```
这样，如果我们在命令行中没有定义 `NAME` 变量，那么 `NAME` 就会被设置为 `xiaoming`。


``` bash
$ make hello
Hello xiaoming

$ make hello NAME=lihua
Hello lihua
```

### 文本替换
在某些情况下，您想从现有字符串创建一个新字符串。

例如，你想从路径中提取目录名或文件名，如

PATH = data/reads/abc.txt
make 语法提供了函数来执行此操作，例如：

DIR = $(dir ${PATH}) 将包含 data/reads/
FNAME = $(notdir ${PATH}) 将包含 abc.txt


``` makefile
# Sets the prefix for commands.
.RECIPEPREFIX = >

# Set a filename
FILE = data/refs/ebola.fa.gz

demo:

# Prints: data/refs/ebola.fa.gz
> @echo ${FILE}

# Prints: data/refs
> @ echo $(dir ${FILE})

# Prints: ebola.fa.gz
> @ echo $(notdir ${FILE})

# Prints: data/refs/human.fa.gz
> @ echo $(subst ebola,human,${FILE})

# Prints: data/refs/ebola.fa
> @echo $(patsubst %.gz,%,${FILE})
```


``` bash
$ make demo
data/refs/ebola.fa.gz
data/refs/
ebola.fa.gz
data/refs/human.fa.gz
data/refs/ebola.fa
```

## 生信流程模块化

虽然我们都希望一个脚本能够完成所有的工作，这样是可实现的，但可能会遇到各种各样的问题，比如脚本过长，不易维护，不易复用等。因此，我们可以将脚本分解为多个模块，每个模块负责一个特定的任务。这样可以提高代码的可读性，可维护性和可复用性。

明白了需求后，按照一个良好的编写习惯，可以很好的将脚本模块化。这样可以提高代码的可读性，可维护性和可复用性。我们接下来编写一个从 SRA 数据库下载 metadata 数据的 make 脚本吧。

### step1: 默认设置

``` makefile
# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory
```

### step2: 打印帮助
我们总是会忘记自己的脚本是干什么的。因此，最好第一个命令是打印帮助。

``` makefile
# Print usage information.
usage:
> @echo "#"
> @echo "# metadata.mk: downloads metadata from SRA"
> @echo "#"
> @echo "# SRA=${SRA}"
> @echo "#"
> @echo "# make run|keys|clean"
> @echo "#"
```

### step3: 定义变量
刚开始可以不需要这一步，等到足够熟练使用后，可以将一些常用的变量定义在这里。

``` makefile
# Sets the default target.
SRA ?= PRJEB31790
```

### step4: 添加代码
这里添加一下下载，提取信息，清理的代码。

``` makefile
run:
> @bio search ${SRA} -H --csv --all > ${SRA}.csv

keys: ${SRA}.csv
> @cat ${SRA}.csv | csvcut -c run_accession,sample_title

clean:
> @rm -f ${SRA}.csv
```

### step5: 运行命令
把上述所有代码放到一个文件中，然后运行 make 命令。

``` bash
# 使用默认变量
make -f metadata.mk run

# 使用自定义变量
make -f metadata.mk run SRA=PRJNA932187

# 提取关键信息
make -f metadata.mk keys
# 也可以直接运行，因为定义了依赖，会默认先运行 run，然后运行 keys

# 清理
make -f metadata.mk clean
```


## 拓展
我经常使用的软件 `bio` :
https://www.bioinfo.help/index.html
这个软件有很多我们常用的生信工具。我建议你使用它！

下载软件后，运行下面这个代码会下载很多的 makefile 脚本，熟悉使用它们，你会发现生信分析变得更加简单。

``` bash
bio code
```



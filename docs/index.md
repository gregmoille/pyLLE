---
layout: index
title: Home
---

# Welcome to pyLLE documentation

_Current version: v3.0.0_


pyLLE is a tool to solve the Lugiato Lefever Equations (LLE) in a fast and easy way. Thanks to a user-friendly front-end in python and an efficient back end in Julia, solving this problem becomes easy and fast.

Using the **open source, free and efficient** computing of Julia, especially the way the FFT is implemented, makes the overall simulation last only a few minutes. On the other hand, python allows for easy scripting, easy display of the figures, and easy saving of the results. 

<i class="fas fa-tachometer-alt"></i> How much faster is pyLLE compared to other implementations? We ran this quick benchmark<sup>[1](#myfootnote1)</sup> to find out. One can see that pure Matlab R2018a is doing not that bad (thanks to a better JIT, beware of previous versions), but in addition to being slower than pyLLE, the major drawback is the proprietary license and low portability of the language. Pure Python is really slow, hence one can see the big addition of doing the intensive calculations in Julia.

| Matlab R2018a <i class="far fa-meh"></i>| Python Only <i class="far fa-sad-cry"></i>|  pyLLE <i class="far fa-thumbs-up"></i> |
|:------:|:-----------:|:-------:|
| 19 min | 45 min |  **11 min**  |

<i class="fas fa-chalkboard-teacher"></i>  For a fairly complete example, please refer to the [example](https://gregmoille.github.io/pyLLE/Example.html) where we investigate the soliton generation in a micro-ring resonator

<i class="far fa-smile-beam"></i> If you are happy with the software, feel free to check out our [paper]() on this package, and [cite us](https://gregmoille.github.io/pyLLE/HowToCite.html) when you publish. 

<i class="far fa-frown-open"></i> If you are not happy with the package (no worry, it happens), please [refer](https://github.com/gregmoille/pyLLE/issues) with any bugs or inquiries to improve the package! Feel free to send us an email if you feel like it (can be found on the [github profile](https://github.com/gregmoille)).

<br>
<br>
<br>
---

<a name="myfootnote1"><sup>1</sup></a> *Simulations made on a MacBook Pro 2017 3.1 GHz Intel Core i5, RAM 16 GB 2133 MHz LPDDR3. System simulated is using the exact same parameters than in the [example](https://gregmoille.github.io/pyLLE/Example.html)*

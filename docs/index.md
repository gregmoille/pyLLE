---
layout: index
title: Home
---

# Welcome to pyLLE documentation

_Current version: v2.1_


pyLLE is a tool to solve the Lugiato Lefever Equations (LLE) in a fast and easy way. Thanks to an user-friendly front-end in python and a efficient back end in Julia, solving this problem becomes easy and fast.

Using the efficient computing of Julia, especially the way fft are implemented makes the overall simulation last only few minutes. On the other hand, python allows a easy scripting, easy display of the figures and easy saving of the results. 

<i class="fas fa-tachometer-alt"></i> How faster pyLLE is compared to other implementation? Well, this quick benchmark<sup>[1](#myfootnote1),</sup> <sup>[2](#myfootnote1)</sup>speaks for itself: 

| Matlab <i class="far fa-dizzy"></i>| Python Only <i class="fas fa-flushed"></i> |  pyLLE <i class="far fa-thumbs-up"></i> |
|:------:|:-----------:|:-------:|
| Matlab | Python Only |  pyLLE  |

For a fairly complete example, please refer to the [example](https://gregmoille.github.io/pyLLE/Example.html) of  which investigate the soliton generation in a micro-ring resonator

<i class="far fa-smile-beam"></i> If you are happy with the software, feel free to check out our [paper]() on this package, and [cite us](https://gregmoille.github.io/pyLLE/HowToCite.html) when you publish. 

<i class="far fa-frown-open"></i> If you are not happy with the package (no worry, it happens), please [refer](https://github.com/gregmoille/pyLLE/issues) any bugs or inquiries to improve the package! 

<br>
<br>
<br>
---

<a name="myfootnote1"><sup>1</sup></a>: *Simulations made on a MacBook Pro 2017 3.1 GHz Intel Core i5, RAM 16 GB 2133 MHz LPDDR3*

<a name="myfootnote2"><sup>2</sup></a>: *System simulated is using the exact same parameter than in the [example](https://gregmoille.github.io/pyLLE/Example.html)*
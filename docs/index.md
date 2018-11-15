---
layout: index
title: Home
---

# Welcome to pyLLE documentation

_Current version: v2.1_


pyLLE is a tool to solve the Lugiato Lefever Equations (LLE) in a fast and easy way. Thanks to an user-friendly front-end in python and a efficient back end in Julia, solving this problem becomes easy and fast.

Using the **open source, free and efficient** computing of Julia, especially the way fft are implemented makes the overall simulation last only few minutes. On the other hand, python allows a easy scripting, easy display of the figures and easy saving of the results. 

<i class="fas fa-tachometer-alt"></i> How faster pyLLE is compared to other implementation? We ran this quick benchmark<sup>[1](#myfootnote1)</sup> to find out. One can see that pure Matlab R2018 is doing well, but the big drawback is the proprietary license and low portability of the langague. Pure Python is really fast, hence one can see the big addition of doing the intensive job in Julia.

| Matlab R2018a <i class="fas fa-meh"></i>| Python Only <i class="far fa-dizzy"></i>| |  pyLLE <i class="far fa-thumbs-up"></i> |
|:------:|:-----------:|:-------:|
|  |  |  **10min 51s**  |

<i class="fas fa-chalkboard-teacher"></i>  For a fairly complete example, please refer to the [example](https://gregmoille.github.io/pyLLE/Example.html) of  which investigate the soliton generation in a micro-ring resonator

<i class="far fa-smile-beam"></i> If you are happy with the software, feel free to check out our [paper]() on this package, and [cite us](https://gregmoille.github.io/pyLLE/HowToCite.html) when you publish. 

<i class="far fa-frown-open"></i> If you are not happy with the package (no worry, it happens), please [refer](https://github.com/gregmoille/pyLLE/issues) any bugs or inquiries to improve the package! Feel free to send us an email if you feel like it (can be found on the [github profile](https://github.com/gregmoille)).

<br>
<br>
<br>
---

<a name="myfootnote1"><sup>1</sup></a> *Simulations made on a MacBook Pro 2017 3.1 GHz Intel Core i5, RAM 16 GB 2133 MHz LPDDR3. System simulated is using the exact same parameters than in the [example](https://gregmoille.github.io/pyLLE/Example.html)*
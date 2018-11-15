---
layout: page
title: pyLLE Modules
sidebar_link: true
---


<div class="section" id="pylle-package">
<h1>pyLLE package<a class="headerlink" href="#pylle-package" title="Permalink to this headline">¶</a></h1>
<div class="section" id="module-pyLLE.analyzedisp">
<span id="pylle-analyzedisp-module"></span><h2>pyLLE.analyzedisp module<a class="headerlink" href="#module-pyLLE.analyzedisp" title="Permalink to this headline">¶</a></h2>
<dl class="class">
<dt id="pyLLE.analyzedisp.AnalyzeDisp">
<em class="property">class </em><code class="descclassname">pyLLE.analyzedisp.</code><code class="descname">AnalyzeDisp</code><span class="sig-paren">(</span><em>**kwargs</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/analyzedisp.html#AnalyzeDisp"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.analyzedisp.AnalyzeDisp" title="Permalink to this definition">¶</a></dt>
<dd><p>Bases: <code class="xref py py-class docutils literal"><span class="pre">object</span></code></p>
<p>Calls to analyze the dispersion of a simulated resonator
Initialization input. Everything is in SI ([]=facultative):
<a href="#id1"><span class="problematic" id="id2">**</span></a>kwargs</p>
<blockquote>
<div><ul class="simple">
<li>f_center &lt;float&gt;: pump frequency</li>
<li>file &lt;str&gt;: .txt file to load</li>
<li>R &lt;float&gt;: radius of the resonator</li>
<li>rM_fit &lt;list&gt;.: lower and upper bonds of mode to fit the dispersion</li>
<li>rM_sim &lt;list&gt;.: lower and upper bonds of the modes to extrapolate for the simulation</li>
<li>f &lt;obj&gt;: matplotlib figure handle</li>
<li>ax &lt;obj&gt;: matplotlib axes handle</li>
<li>label &lt;list&gt;: list of the string for each plot to be labelled</li>
<li>plottype &lt;str&gt;: define the type of plot &#8216;all&#8217; [defaults], &#8216;sim&#8217;, &#8216;fit&#8217;, &#8216;fem&#8217;, &#8216;ind&#8217;</li>
</ul>
</div></blockquote>
<dl class="method">
<dt id="pyLLE.analyzedisp.AnalyzeDisp.GetDint">
<code class="descname">GetDint</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/analyzedisp.html#AnalyzeDisp.GetDint"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.analyzedisp.AnalyzeDisp.GetDint" title="Permalink to this definition">¶</a></dt>
<dd><p>Retrieve the dispersion of a resonator based on the frequency of 
resonance and azimuthal mode order. The data are fit using a cubic spline method</p>
<p><strong>Output</strong></p>
<blockquote>
<div><ul class="simple">
<li>self.PrM_fit: scipy.interpolate object which fitted the data</li>
<li>self.Dint_fit: fitted integrated dispersion for the simulated mode</li>
<li>self.neff_pmp: effective index at the pump frequency</li>
<li>self.ng_pmp: group index at the pump frequency</li>
</ul>
</div></blockquote>
</dd></dl>

</dd></dl>

</div>
<div class="section" id="module-pyLLE.llesolver">
<span id="pylle-llesolver-module"></span><h2>pyLLE.llesolver module<a class="headerlink" href="#module-pyLLE.llesolver" title="Permalink to this headline">¶</a></h2>
<dl class="class">
<dt id="pyLLE.llesolver.LLEsovler">
<em class="property">class </em><code class="descclassname">pyLLE.llesolver.</code><code class="descname">LLEsovler</code><span class="sig-paren">(</span><em>**kwargs</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler" title="Permalink to this definition">¶</a></dt>
<dd><p>Bases: <code class="xref py py-class docutils literal"><span class="pre">object</span></code></p>
<p>Class to solve the Lugiato Lefever Equation
Initialization input ([]=facultative):</p>
<p><strong>res &lt;dict&gt;</strong></p>
<blockquote>
<div><ul class="simple">
<li>Qi &lt;float&gt;: intrinsic Q of the resonator</li>
<li>Qc &lt;float&gt;: coupling Q of the resonator</li>
<li>R &lt;float&gt;: ring radius</li>
<li>gamma &lt;float&gt;: Non linear index of the material</li>
<li>dispfile &lt;str&gt; : str pointing to a .txt file where the azimuthal mode orders and corresponding resonances are saved and separated by a coma</li>
</ul>
</div></blockquote>
<p><strong>sim &lt;dict&gt;</strong></p>
<blockquote>
<div><ul class="simple">
<li>Tscan &lt;float&gt;: length of the simulation (in unit of round trip)</li>
<li>mu_fit &lt;list&gt;: number of mode to fit</li>
<li>mu_sim &lt;list&gt;: number of mode to simulate</li>
<li>domega_init &lt;float&gt;: initial detuning of the pump</li>
<li>domega_end &lt;float&gt;: final detuning of the pump</li>
</ul>
</div></blockquote>
<p><strong>debug &lt;bool&gt;</strong>: Save a trace in a log file of the different actions pyLLE perform (default = True)</p>
<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.Analyze">
<code class="descname">Analyze</code><span class="sig-paren">(</span><em>plot=False</em>, <em>f=None</em>, <em>ax=None</em>, <em>label=None</em>, <em>plottype='all'</em>, <em>zero_lines=True</em>, <em>mu_sim=None</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.Analyze"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.Analyze" title="Permalink to this definition">¶</a></dt>
<dd><p>Call pyLLE.analyzedisp.AnalyzeDisp to get the dispersion of the resonator</p>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.PlotCombPower">
<code class="descname">PlotCombPower</code><span class="sig-paren">(</span><em>do_matplotlib=False</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.PlotCombPower"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.PlotCombPower" title="Permalink to this definition">¶</a></dt>
<dd><p>Plot a figure with 3 subplots.</p>
<ul class="simple">
<li>Top subplot = map of the spectra for the steps taken by the LLE (step sub-sampled to be 1000)</li>
<li>middle subplot = temporal map of the intensity inside the resonator for the steps of the LLE</li>
<li>bottom subplot = normalized comb power</li>
</ul>
<p><strong>Output</strong></p>
<blockquote>
<div><ul class="simple">
<li><dl class="first docutils">
<dt>if jupyter notebook:</dt>
<dd><ul class="first last">
<li>fig:  handle of figure and axes of the plotly figure displayed</li>
</ul>
</dd>
</dl>
</li>
<li><dl class="first docutils">
<dt>if ipython or script:</dt>
<dd><ul class="first last">
<li>f, ax:  handle of figure and axes of the matplotlib figure displayed</li>
</ul>
</dd>
</dl>
</li>
</ul>
</div></blockquote>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.PlotCombSpectra">
<code class="descname">PlotCombSpectra</code><span class="sig-paren">(</span><em>ind</em>, <em>f=None</em>, <em>ax=None</em>, <em>label=None</em>, <em>pwr='both'</em>, <em>do_matplotlib=False</em>, <em>plot=True</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.PlotCombSpectra"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.PlotCombSpectra" title="Permalink to this definition">¶</a></dt>
<dd><p>Plot the spectra for a given index in the 1000 sub-sampled LLE steps</p>
<p>It save the corresponding spectra in the attribute <strong>self.spectra</strong></p>
<p><strong>Input</strong></p>
<blockquote>
<div><ul class="simple">
<li>ind &lt;ind&gt;: index in the LLE step to plot the spectra</li>
<li>f &lt;obj&gt;:  matplotlib figure handle (if None, new figure)</li>
<li>ax &lt;obj&gt;: matplotlib axe handle</li>
<li>label &lt;str&gt;: label for the legend</li>
<li>pwr &lt;str&gt;: &#8216;both&#8217;, &#8216;ring&#8217;, &#8216;wg&#8217; depending on the spectra wanted (inside the ring, the waveguide or both)</li>
</ul>
</div></blockquote>
<p><strong>Output</strong></p>
<blockquote>
<div><blockquote>
<div><ul class="simple">
<li><dl class="first docutils">
<dt>if jupyter notebook:</dt>
<dd><ul class="first last">
<li>fig:  handle of figure and axes of the plotly figure displayed</li>
</ul>
</dd>
</dl>
</li>
</ul>
</div></blockquote>
<ul class="simple">
<li><dl class="first docutils">
<dt>if ipython or script:</dt>
<dd><ul class="first last">
<li>f, ax:  handle of figure and axes of the matplotlib figure displayed</li>
</ul>
</dd>
</dl>
</li>
</ul>
</div></blockquote>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.PlotSolitonTime">
<code class="descname">PlotSolitonTime</code><span class="sig-paren">(</span><em>ind</em>, <em>f=None</em>, <em>ax=None</em>, <em>label=None</em>, <em>do_matplotlib=False</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.PlotSolitonTime"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.PlotSolitonTime" title="Permalink to this definition">¶</a></dt>
<dd><p>Plot the temporal profile for a given index in the 1000 sub-sampled LLE step</p>
<p>It save the corresponding time profile in the attribute <strong>self.fasttime</strong></p>
<p><strong>Input</strong></p>
<blockquote>
<div><ul class="simple">
<li>ind &lt;ind&gt;: index in the LLE step to plot the spectra</li>
<li>f &lt;obj&gt;:  matplotlib figure handle (if None, new figure)</li>
<li>ax &lt;obj&gt;: matplotlib axe handle</li>
<li>label &lt;str&gt;: label for the legend</li>
</ul>
</div></blockquote>
<p><strong>Output</strong></p>
<blockquote>
<div><blockquote>
<div><ul class="simple">
<li><dl class="first docutils">
<dt>if jupyter notebook:</dt>
<dd><ul class="first last">
<li>fig:  handle of figure and axes of the plotly figure displayed</li>
</ul>
</dd>
</dl>
</li>
</ul>
</div></blockquote>
<ul class="simple">
<li><dl class="first docutils">
<dt>if ipython or script:</dt>
<dd><ul class="first last">
<li>f, ax:  handle of figure and axes of the matplotlib figure displayed</li>
</ul>
</dd>
</dl>
</li>
</ul>
</div></blockquote>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.RetrieveData">
<code class="descname">RetrieveData</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.RetrieveData"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.RetrieveData" title="Permalink to this definition">¶</a></dt>
<dd><p>Load the output hdf5 saved by julia and transform it in a user-friendly dictionary to be more pythonistic</p>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.SavePlots2File">
<code class="descname">SavePlots2File</code><span class="sig-paren">(</span><em>basename='./'</em>, <em>format='pdf'</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.SavePlots2File"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.SavePlots2File" title="Permalink to this definition">¶</a></dt>
<dd><p>Save the different figures in the corresponding format with matplotlib and the previous class Latexify. It will only save the figure of the type of plot displayed earlier
The name of the file are automatically created depending of the plot</p>
<p><strong>Input</strong></p>
<blockquote>
<div><ul class="simple">
<li>basename &lt;str&gt;: basename of the file. It can include a path too</li>
<li>format &lt;str&gt;: format to save the figure. Has to be compatible with matplotlib</li>
</ul>
</div></blockquote>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.SaveResults">
<code class="descname">SaveResults</code><span class="sig-paren">(</span><em>fname</em>, <em>path='./'</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.SaveResults"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.SaveResults" title="Permalink to this definition">¶</a></dt>
<dd><p>Save the whole class with pickle to be able to easilly call it back or retrieve the results after saving</p>
<p><strong>Input</strong></p>
<blockquote>
<div><ul class="simple">
<li>fname &lt;str&gt;: name to save. The &#8216;.pkl&#8217; extension will be added</li>
<li>path &lt;str&gt;: path to save the results (defaults &#8216;./&#8217;)</li>
</ul>
</div></blockquote>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.Setup">
<code class="descname">Setup</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.Setup"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.Setup" title="Permalink to this definition">¶</a></dt>
<dd><p>Setup the simulation for the Julia back-end. 
Save the two main dictionary self.sim and self.res into a readable hdf5 file for Julia in the temporary location define by the os</p>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.SolveSteadySteate">
<code class="descname">SolveSteadySteate</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.SolveSteadySteate"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.SolveSteadySteate" title="Permalink to this definition">¶</a></dt>
<dd><p>Newton Method to find the root of the steady state equation</p>
</dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.LLEsovler.SolveTemporal">
<code class="descname">SolveTemporal</code><span class="sig-paren">(</span><em>tol=0.001</em>, <em>maxiter=6</em>, <em>step_factor=0.1</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#LLEsovler.SolveTemporal"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.SolveTemporal" title="Permalink to this definition">¶</a></dt>
<dd><p>Call Julia to solve the LLE. 
Julia save in a temp file the progress and this method display it through a progress bar.
The subprocess.Popen is saved in the attribute self.JuliaSolver</p>
</dd></dl>

<dl class="attribute">
<dt id="pyLLE.llesolver.LLEsovler.hbar">
<code class="descname">hbar</code><em class="property"> = 1.0558338924716337e-34</em><a class="headerlink" href="#pyLLE.llesolver.LLEsovler.hbar" title="Permalink to this definition">¶</a></dt>
<dd></dd></dl>

</dd></dl>

<dl class="class">
<dt id="pyLLE.llesolver.Latexify">
<em class="property">class </em><code class="descclassname">pyLLE.llesolver.</code><code class="descname">Latexify</code><span class="sig-paren">(</span><em>**kwarg</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#Latexify"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.Latexify" title="Permalink to this definition">¶</a></dt>
<dd><p>Bases: <code class="xref py py-class docutils literal"><span class="pre">object</span></code></p>
<p>Class that handle saving the figures in a nice way compatible with 
the column/page size of different latex template</p>
<dl class="docutils">
<dt>input [] = optional: </dt>
<dd><ul class="first last simple">
<li>figname = name to save the figure (without extension)</li>
<li></li>
<li>fig = matplotlib handle to the figure</li>
<li>[fig_width]: default = 1column</li>
<li>[frmt]: default = pdf</li>
<li>[fig_height] : default = 6.5</li>
<li>[font_size] : default = 8</li>
</ul>
</dd>
</dl>
<dl class="method">
<dt id="pyLLE.llesolver.Latexify.GoBackToNormal">
<code class="descname">GoBackToNormal</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#Latexify.GoBackToNormal"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.Latexify.GoBackToNormal" title="Permalink to this definition">¶</a></dt>
<dd></dd></dl>

<dl class="method">
<dt id="pyLLE.llesolver.Latexify.SavePlot">
<code class="descname">SavePlot</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#Latexify.SavePlot"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.Latexify.SavePlot" title="Permalink to this definition">¶</a></dt>
<dd></dd></dl>

</dd></dl>

<dl class="class">
<dt id="pyLLE.llesolver.MyLogger">
<em class="property">class </em><code class="descclassname">pyLLE.llesolver.</code><code class="descname">MyLogger</code><span class="sig-paren">(</span><em>fname</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#MyLogger"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.MyLogger" title="Permalink to this definition">¶</a></dt>
<dd><p>Bases: <code class="xref py py-class docutils literal"><span class="pre">object</span></code></p>
<p>Custom made logger as the logger default package cannot be pickled</p>
<dl class="method">
<dt id="pyLLE.llesolver.MyLogger.info">
<code class="descname">info</code><span class="sig-paren">(</span><em>method</em>, <em>message</em><span class="sig-paren">)</span><a class="reference internal" href="_modules/pyLLE/llesolver.html#MyLogger.info"><span class="viewcode-link">[source]</span></a><a class="headerlink" href="#pyLLE.llesolver.MyLogger.info" title="Permalink to this definition">¶</a></dt>
<dd></dd></dl>

</dd></dl>

</div>
</div>


           </div>
           
          </div>
          <footer>
  
    <div class="rst-footer-buttons" role="navigation" aria-label="footer navigation">
      
      
        <a href="team.html" class="btn btn-neutral" title="pyLLE Team" accesskey="p" rel="prev"><span class="fa fa-arrow-circle-left"></span> Previous</a>
      
    </div>
  

  <hr/>

  <div role="contentinfo">
    <p>
        &copy; Copyright 2018, Gregory Moille

    </p>
  </div>
  Built with <a href="http://sphinx-doc.org/">Sphinx</a> using a <a href="https://github.com/rtfd/sphinx_rtd_theme">theme</a> provided by <a href="https://readthedocs.org">Read the Docs</a>. 

</footer>

        </div>
      </div>

    </section>

  </div>
  


  

    
    
      <script type="text/javascript">
          var DOCUMENTATION_OPTIONS = {
              URL_ROOT:'./',
              VERSION:'2.0',
              LANGUAGE:'None',
              COLLAPSE_INDEX:false,
              FILE_SUFFIX:'.html',
              HAS_SOURCE:  true,
              SOURCELINK_SUFFIX: '.txt'
          };
      </script>
        <script type="text/javascript" src="_static/jquery.js"></script>
        <script type="text/javascript" src="_static/underscore.js"></script>
        <script type="text/javascript" src="_static/doctools.js"></script>
        <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
    

  

  <script type="text/javascript" src="_static/js/theme.js"></script>

  <script type="text/javascript">
      jQuery(function () {
          SphinxRtdTheme.Navigation.enable(true);
      });
  </script> 

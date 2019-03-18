---
author: Remco Bouckaert
level: Intermediate
title: Model selection
subtitle: with nested sampling
beastversion: 2.5.2
nsversion: 1.0.3
---


# Background

BEAST provides a bewildering number of models. Bayesians have two techniques to deal with this problem: model averaging and model selection. With **model averaging**, the MCMC algorithms jumps between the different models, and models more appropriate for the data will be sampled more often than unsuitable models (see for example the [substitution model averaging tutorial](https://taming-the-beast.org/tutorials/Substitution-model-averaging/)). **Model selection** on the other hand just picks one model and uses only that model during the MCMC run. This tutorial gives some guidelines on how to select the model that is most appropriate for your analysis. 

Bayesian model selection is based on estimating the marginal likelihood: the term forming the denominator in Bayes formula. This is generally a computationally intensive task and  there are several ways to estimate them. Here, we concentrate on nested sampling as a way to estimate the marginal likelihood as well as the uncertainty in that estimate.

Say, we have two models, M1 and M2, and estimates of the (log) marginal likelihood, ML1 and ML2, then we can calculate the Bayes factor, which is the fraction BF=ML1/ML2 (or in log space, the difference log(BF) = log(ML1)-log(ML2)). If BF is larger than 1, model M1 is favoured, and otherwise M2 is favoured. How much it is favoured can be found in the following table ({% cite kass1995bayes --file Tutorial-Template/master-refs.bib %})

<img style="width:80.0%;" src="figures/BFs.png" alt="">

Note that sometimes a factor 2 is used for multiplying BFs, so when comparing BFs from different publications, be aware which definition that was used.


----

# Programs used in this tutorial 

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014,bouckaert2018beast --file Tutorial-Template/master-refs.bib %}. This tutorial uses the BEAST2 version 2.5.2.

### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### Tracer

[Tracer](http://tree.bio.ed.ac.uk/software/tracer) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v1.7.0.

----

# Practical: Selecting a clock model

We will analyse a set of hepatitis B virus (HBV) sequences samples through time and concentrate on selecting a clock model. The most popular clock models are the strict clock model and uncorrelated relaxed clock with log normal distributed rates (UCLN) model.

## Setting up the analyses

First thing to do is set up the two analyses in BEAUti, and run them in order to make sure there are differences in the analyses. In BEAUti:

* start a new analysis using the Standard template
* import HBV.nex (menu File/Import alignment)
* In the tip-dates panel, select "tip dates", click "Auto configure" and select the "split on character" option, taking group 2 (see [Fig 1](#fig:auto-config)).
* in the site model panel, select HKY as substitution model and leave the rest as is.
* in the clock model panel, set the clock rate to 2e-5. Though usually, we want to estimate the rate, to speed things up for the tutorial, we fix the clock rate at that number as follows:
    * uncheck menu "Mode/Automatic set clock rate". Now the estimate entry should not be grayed out any more.
    * uncheck the "estimate" box next to the clock rate entry.
* in the priors panel, select "Coalescent Constant Population" as tree prior.
* also in the priors panel, change to popSize prior to Gamma with alpha = 0.01, beta = 100 ([Fig 2](#fig:prior))
* in the MCMC panel, change the Chain Length to 1 million.
* you can rename the file for trace log and tree file to include "Strict" to distinguish them for the relaxed clock ones.
* save the file as HBVStrict.xml
* run the analysis with BEAST

<figure>
	<a id="fig:auto-config"></a>
	<img style="width:80.0%;" src="figures/BEAUti-configure-tip-dates.png" alt="">
	<figcaption>Figure 1: Configuring tip dates in BEAUti</figcaption>
</figure>
<br>

<figure>
	<a id="fig:priors"></a>
	<img style="width:80.0%;" src="figures/BEAUti-priors.png" alt="">
	<figcaption>Figure 2: Priors panel for strict clock analysis in BEAUti</figcaption>
</figure>
<br>

While you are waiting for BEAST to finish, it is time to set up the relaxed clock analysis. This is now straightforward:
* in the clock model panel, change "Strict clock" to "Relaxed Clock Log Normal".
* set the clock rate to 2e-5, and uncheck the "estimate" box.
* in the MCMC panel, replace "strict" in the file names for trace and tree log to "tree".
* save file as HBVUCLN.xml
* run the analysis in BEAST

Once the analyses have run, open the log file in Tracer and compare estimates and see whether the analyses substantially differ. You can also compare the trees in DensiTree.

> Are there any statistics that are very different?
> Do tree heights differ substantially? 
> Which analysis is preferable and why? 
<!-- depends on the question you want to answer: if tree height is of interest, strict clock is preferred, since it reduces the uncertainty. If kappa is of interest, things are not that different -->

If there are no substantial differences between the analysis for the question you are interested in, you do not have to commit to one model or another, and you can claim that the results are robust under different models. However, if there are significant differences, you may want to do a formal test to see which model is preferred over other models. In a Bayesian context, in practice this comes down to estimating the marginal likelihood, and calculating Bayes factors: the ratios of marginal likelihoods. Nested sampling {% cite russel2018model --file Tutorial-Template/master-refs.bib %} is one way to estimate marginal likelihoods. 

## Installing the NS Package

To use nested sampling, first have to install the NS (version {{ page.nsversion }} or above) package.

> Open BEAUti and navigate to **File > Manage Packages**. Select MODEL_SELECTION and then click **Install/Upgrade** ([Fig 3](#fig:install)). Then **_restart BEAUti_** to load the package.
>  

<figure>
	<a id="fig:install"></a>
	<img style="width:80.0%;" src="figures/install_NS.png" alt="">
	<figcaption>Figure 3: Installing NS in the Manage Packages window in BEAUti</figcaption>
</figure>
<br>


-------

# Tutorial style guide

## Text styling

This is how to write _italic text_.

This is how to write **bold text**.

This is how to write **_bold and italic text_**.

Do text superscripts like this 7^th, x^2y or  x^(2y + 3z).


## Lists

### Unnumbered lists

- Lorem ipsum dolor sit amet, consectetur adipiscing elit.
- Integer pharetra arcu ut nisl mollis ultricies.
	- Fusce nec tortor at enim cursus dictum.
	- Phasellus nec urna quis velit eleifend convallis sodales nec augue.
- In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
- Nam vitae turpis eu lacus imperdiet mollis id at augue.
- Sed sed turpis ac dolor mollis accumsan.


### Numbered lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	1. Fusce nec tortor at enim cursus dictum.
	2. Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.

### Mixed lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	* Fusce nec tortor at enim cursus dictum.
	* Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.


## Figures


<figure>
	<a id="fig:example1"></a>
	<img style="width:25%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 1: This figure is 25% of the page width.</figcaption>
</figure>


<figure>
	<a id="fig:example2"></a>
	<img style="width:10%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 2: This figure is only 10% of the page width.</figcaption>
</figure>



# Code

A bit of inline monospaced font can be made `like this`. Larger code blocks can be made by using the code environment:

Java:

```java
public class HelloWorld {

    public static void main(String[] args) {
        // Prints "Hello, World" to the terminal window.
        System.out.println("Hello, World");
    }

}
```

XML:

```xml
	<BirthDeathSkylineModel spec="BirthDeathSkylineModel" id="birthDeath" tree="@tree" contemp="true">
	      <parameter name="origin" id="origin" value ="100" lower="0."/>    
	      <parameter name="R0" id="R0" value="2" lower="0." dimension ="10"/>
	      <parameter name="becomeUninfectiousRate" id="becomeUninfectiousRate" value="1" lower="0." dimension ="10"/>
	      <parameter name="samplingProportion" id="samplingProportion" value="0."/>
	      <parameter name="rho" id="rho" value="1e-6" lower="0." upper="1."/>
	</BirthDeathSkylineModel>
```

R:

```R
	> myString <- "Hello, World!"
	> print (myString)
	[1] "Hello, World!"
```

# Equations

Inline equations: {% eqinline \dot{x} = \sigma(y-x) %}

Displayed equations: 
{% eq \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right) %}



## Instruction boxes

Use block-quotes for step-by-step instruction that the user should perform (this will produce a framed box on the website):

> The data we have is not the data we want, and the data we need is not the data we have.
> 
> We can input **any** formatted text in here:
>
> - Even
> - Lists
>
> or equations:
>
> {% eq (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right) %}






# Hyperlinks

Add links to figures like this: 

- [Figure 1](#fig:example1) is 25% of the page width.
- [Figure 2](#fig:example2) is 10% of the page width. 

Add links to external URLs like [this](http://www.google.com). 

Links to equations or different sections within the same document are a little buggy.


----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Tutorial-Template/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Tutorial-Template/master-refs.bib %}


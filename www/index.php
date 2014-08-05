<?php
include 'includes/header.php';
?>

<body>  

<section id=sec_summary>
	<h1 onclick="layout_showunshow('div_summary',this,'inline-block');">Summary &#x25BE;</h1>
	<div id=div_summary style='display:inline-block;'>

		<i>ProbeAnnotator</i> is a package that contains simple and efficient functions for the annotation of genomic platforms (e.g. Illumina, Affymetrix) using bowtie alignment.<br> This tool as been developed within the <i>Bioinformatics Core Facility, Swiss Institute of Bioinformatics</i>. More information about the team on our <a href="http://bcf.isb-sib.ch/">webpage</a>.<br> 
		<p> The project summary page on <strong>r-forge</strong> can be found <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/">here</a>.</p>

	</div> 
</section> 


<section id=sec_tutorial>
	<h1 onclick="layout_showunshow('div_tutorial',this,'inline-block');">Tutorials &#x25BE;</h1>
	<div id=div_tutorial style='display:inline-block;'>
		<table>
		<thead>
		<tr><th>Name</th><th>Platform</th><th>Version</th><th>Files</th></tr>
		</thead>
		<tbody>
		</tbody>
		</table>
	</div> 
</section> 



<section id=sec_download>
	<h1 onclick="layout_showunshow('div_download',this,'inline-block');">Download Package &#x25BE;</h1>
	<div id=div_download style='display:inline-block;'>
		<table>
		<thead>
		<tr><th>System</th><th>Version</th><th>File</th></tr>
		</thead>
		<tbody>
		<tr><td>Linux</td><td>1.0.1</td><td><a href=ProbeAnnotator_1.0.1.tar.gz>download</td></tr>
		</tbody>
		</table>
	</div> 
</section> 



<section id=sec_contact>
	<h1 onclick="layout_showunshow('footer',this,'block');">Contact &#x25BE;</h1>
	<div id=footer style='display:block;'>
		<p id="about">If any information is needed, contact the (small) developer team.</p>
		<div id="footerlist"> 
			<h3>Contact</h3>
			<a href="mailto:alexandre.thiery@unil.ch">Via Email</a>
		</div> 
		<div id=team> 
			<h3>Team</h3>
			<p id="footercopy">Alexandre Thiéry, Swiss Institute of Bioinformatics</p>
		</div> 
	</div> 
</section> 
	

</html>


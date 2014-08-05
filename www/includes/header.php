<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html>
<html>
<head>  
	<!--<script language="javascript" type="text/javascript" src="../js/jquery.js"></script>--> 
	<meta charset="utf-8">  
	<title>Welcome to ProbeAnnotator</title>  
	

    <meta name="description" content="ProbeAnnotator is a package that contains simple and efficient functions for the annotation of genomic platforms.">
	<link rel="stylesheet" href="styles/custom.css">    
	<link rel="stylesheet" href="styles/tables.css">    
</head>  
<script>
	function layout_showunshow(elemToChange_id,elemEven,display) {
		var str_temp1 = elemEven.innerHTML;

		var res = document.getElementById(elemToChange_id);
		if(res.style.display == display) {
			res.style.display = 'none';
			str_temp1 = str_temp1.replace("▾","◂"); 
		} else {
			res.style.display = display;
			str_temp1 = str_temp1.replace("◂","▾");
		}
		elemEven.innerHTML = str_temp1;
	}
</script>

<header>   
	<p><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /></p>
	<h1>Welcome to <i>ProbeAnnotator</i> project's page.</h1>
</header>



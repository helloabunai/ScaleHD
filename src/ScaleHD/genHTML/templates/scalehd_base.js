/////////////
// GLOBALS //
/////////////
var scaleHDVersion = $('body').data('shd_version')
var instanceLabel = $('body').data('instance_label')
var allParents = document.getElementsByClassName('sequence_sample_link')

///////////////
// FUNCTIONS //
///////////////
$(document).ready(function(){
	// Only for document load events
	// hide #help on load, so that #welcome is the landing page
	$('#help').hide();
	CAGSummaryChart(instanceLabel);
	CCGSummaryChart(instanceLabel);
});

// Function for clicking on a sequence sample
function showSequence(identifier){
	$('div.sequence').hide();
	$('#welcome').hide();
	var seqHead = $('div.sequence[data-sequenceid='+ identifier +']')
	seqHead.show();
}

// Function to change HTML header text if sample/substage link pressed
function changeHeader(replaceRight){
	$('#rightHeader').html(replaceRight)
}

// Function to render lineGraphs within samples when called upon
function lineGraph(identifier, suffix)
{
	var target_id = identifier + suffix;
	var target_element = document.getElementById(target_id)
	var ctx = target_element.getContext('2d');

	console.log(target_element.attributes)

	var chart = new Chart(ctx, {
	  type: 'line',
	  data: {
	    labels: [1500,1600,1700,1750,1800,1850,1900,1950,1999,2050],
	    datasets: [{
	        data: [86,114,106,106,107,111,133,221,783,2478],
	        label: "Africa",
	        borderColor: "#3e95cd",
	        fill: false
	      }, {
	        data: [282,350,411,502,635,809,947,1402,3700,5267],
	        label: "Asia",
	        borderColor: "#8e5ea2",
	        fill: false
	      }, {
	        data: [168,170,178,190,203,276,408,547,675,734],
	        label: "Europe",
	        borderColor: "#3cba9f",
	        fill: false
	      }, {
	        data: [40,20,10,16,24,38,74,167,508,784],
	        label: "Latin America",
	        borderColor: "#e8c3b9",
	        fill: false
	      }, {
	        data: [6,3,2,2,7,26,82,172,312,433],
	        label: "North America",
	        borderColor: "#c45850",
	        fill: false
	      }
	    ]
	  },
	  options: {
	    title: {
	      display: true,
	      text: 'World population per region (in millions)'
	    }
	  }
	});
}

// Function to clear sample links back to white when called
// Clear all sample links of any green colouring
// I.e. un-do changes from the PARENT COLOURING within sublink click event
function clearSampleColour(passedIdentifier){
	for (var i=0; i < allParents.length; i++){
		// If currentIdent of loop matches the clicked identifier
		// clear ALL OTHER idents of any hard-coded HTML styling attributes
		// allows default css style rules (i.e. hover green) to be enforced
		var currIdent = allParents[i].getAttribute('data-sequenceid');
		if (currIdent != passedIdentifier){
			allParents[i].removeAttribute('style');
		}
	}
}

// Function to scroll list to top
function scrollToTop(){
	$('#sideNav').animate({scrollTop: 0}, 'fast')
}

///////////////////////
// LINK CLICK EVENTS //
///////////////////////

// When sequence sample link pressed, call function to show sequence(ID)
$('.sequence_sample_link').click(function(event){
	$('#help').hide();
	$('#welcome').hide();
	event.preventDefault();
	var identifier = $(event.target).data('sequenceid');
	clearSampleColour(identifier)
	showSequence(identifier);

	// render FastQC PBSQ graph
	lineGraph(identifier, '_FQC_PBSQ');
	// render FastQC PBNC graph
	lineGraph(identifier, '_FQC_PBNC');
	// render FastQC seqlen graph
	lineGraph(identifier, '_FQC_SEQLENDIST');

	// Update header to be "SAMPLE | LANDING_STRING"
	rightCandidate = "\
	<p class=\"alignleft\">" + identifier + " | ScaleHD results" + "</p>\n\
	<p class=\"alignright\"><a class=\"contact\" href=\"mailto:alastair.maxwell@glasgow.ac.uk\"><strong>@</strong></a>  |  <a class=\"help\" href=\"#\" data-title=\"Help\"><strong>?</strong></a></p>\n\
	<div style=\"clear: both;\"></div>\n\
	"

	changeHeader(rightCandidate)
	$('html,body').animate({scrollTop: 0},'fast');
});

// When the #home div button is pressed
$('.home').click(function(event){
	clearSampleColour('None')
	$('div.sequence').hide();
	$('#help').hide();
	$('#welcome').show();

	rightCandidate = "\
	<p class=\"alignleft\">ScaleHD - Automated Huntington Disease genotyping</p>\n\
	<p class=\"alignright\"><a class=\"contact\" href=\"mailto:alastair.maxwell@glasgow.ac.uk\"><strong>@</strong></a>  |  <a class=\"help\" href=\"#\" data-title=\"Help\"><strong>?</strong></a></p>\n\
	<div style=\"clear: both;\"></div>\n\
	"
	changeHeader(rightCandidate)
});

// When the #help div button is pressed
// UPDATED -- bind the click event to ANY element with .help -- regardless of when it was made.
$(document).on('click', '.help', function(e){
	$('div.sequence').hide();
	$('#welcome').hide();
	$('#help').show();
	clearSampleColour('None')

	rightCandidate = "\
	<p class=\"alignleft\">ScaleHD - Help!</p>\n\
	<p class=\"alignright\"><a class=\"contact\" href=\"mailto:alastair.maxwell@glasgow.ac.uk\"><strong>@</strong></a>  |  <a class=\"help\" href=\"#\" data-title=\"Help\"><strong>?</strong></a></p>\n\
	<div style=\"clear: both;\"></div>\n\
	"
	changeHeader(rightCandidate)
});

// Sequence sublink (seqqc/aln/gtype) pressed
// Show the parent sequence sample div and scroll to the sublink area
$('.sequence_sample_sublink').click(function(event)
{
	//Function-wide variables
	var identifier = $(event.target).data('sequenceid');
	var title = $(event.target).data('title');
	var longform = "";

	// Clear all sample links of any green colouring
	// I.e. un-do changes from the PARENT COLOURING below
	clearSampleColour(identifier)

	// Create longform string for rightHeader
	if (title === 'seqqc'){
		longform = 'Sequence quality control results'
	}
	else if (title === 'seqalign'){
		longform = 'Sequence alignment results'
	}
	else if (title === 'genotype'){
		longform = 'Allele genotype results'
	}

	//Actually show the sequence on the rightPanel
	showSequence(identifier);
	// render FastQC PBSQ graph
	lineGraph(identifier, '_FQC_PBSQ');
	// render FastQC PBNC graph
	lineGraph(identifier, '_FQC_PBNC');
	// render FastQC seqlen graph
	lineGraph(identifier, '_FQC_SEQLENDIST');

	// PARENT COLOURING
	// If the user directly selects a sub-stage sublink..
	// ..colour the parent (sample name) link also
	targetParent = $(event.target).parent().parent().parent()
	if (targetParent.has('.a.sequence_sample_link')){
		var targetChild = targetParent.children()[0];
		targetChild.style.setProperty('color','#009419');

	}

	// Modifying rightHeader text to be relevant to user selected stage
	rightCandidate = "\
	<p class=\"alignleft\">" + identifier + "<p id=\"seqStage\">   #" + title + "</p><p id=\"longForm\"> | " + longform + "</p>\n\
	<p class=\"alignright\"><a class=\"contact\" href=\"mailto:alastair.maxwell@glasgow.ac.uk\"><strong>@</strong></a>  |  <a class=\"help\" href=\"#\" data-title=\"Help\"><strong>?</strong></a></p>\n\
	<div style=\"clear: both;\"></div>\n\
	"
	changeHeader(rightCandidate)

	// Animate scrolling to position of div
	var header = $('h3[data-sequenceid=' + identifier + '][data-title=' + title + ']');
	$('html,body').animate({scrollTop: header.offset().top - 50},'slow');
});

////////////////////
// Summary Charts //
////////////////////

// CAG Summary
{CAG_FUNCTION}

// CCG Summary
{CCG_FUNCTION}

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
	barGraph(instanceLabel, 'CAGSummaryChart');
	barGraph(instanceLabel, 'CCGSummaryChart');
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

// Function to render barGraphs within summaries and samples when called
function barGraph(identifier, triplet)
{
	var target_id = identifier + triplet;
	var target_element = document.getElementById(target_id)
	var ctx = target_element.getContext('2d');
	var gradient = ctx.createLinearGradient(0,0,0,600);
	gradient.addColorStop(0, '#009419');
	gradient.addColorStop(1, '#001504');

	// get data from JSON attribute
	var data_title = target_element.getAttribute('data-title')
	var data_labels = target_element.getAttribute('data-labels')
	var data_values = target_element.getAttribute('data-values')
	// replace single quotes from python strings to double quotes (JSON A BITCH)
	data_labels = data_labels.replace(/'/g, '"')
	data_values = data_values.replace(/'/g, '"')

	var chart = new Chart(ctx,
	{
		type: 'bar',
		data:
		{
				labels: JSON.parse(data_labels),
				datasets:
				[{
						label: '# of alleles present',
						data: JSON.parse(data_values),
						backgroundColor: gradient,
						hoverBackgroundColor: gradient,
						hoverBorderWidth: 2,
						hoverBorderColor: '#001504'
				}]
		},
		options:
		{
			// Title colour/text
			title:
			{
					display: true,
					fontColor: '#000000',
					fontSize: 15,
					fontFamily: 'Product Sans',
					text: data_title
			},
			// Prevent legend from rendering
			legend:
			{
				display: false
			},
			// X/Y Axes settings
			scales:
			{
				//Y Axis
				yAxes:
				[{
					gridlines:
					{
						display: true
					},
					ticks:
					{
						beginAtZero:false,
						precision:0
					}
				}],
				//X Axis
				xAxes:
				[{
					gridlines:
					{
						display: false,
					},
					ticks:
					{
						autoSkip: true,
						beginAtZero: false,
						stepSize: 20,
						maxTicksLimit: 11
					}
				}]
			}
		 }
	});
}

// Function to render lineGraphs within samples when called upon
function lineGraph(identifier, suffix)
{
	var target_id = identifier + suffix;
	var target_element = document.getElementById(target_id)
	var ctx = target_element.getContext('2d');

	//data from chart attributes since i'm too fuckin dumb to pass it
	//any other way!! fuck web servers!! apparently!!
	var data_title = target_element.getAttribute("data-title")
	var data_labels = target_element.getAttribute("data-labels")
	var data_values = target_element.getAttribute("data-values")
	var data_descr = target_element.getAttribute("data-descr")
	var data_xaxis = target_element.getAttribute("data-xaxis")
	var data_yaxis = target_element.getAttribute("data-yaxis")
	// replace single quotes from python strings to double quotes (JSON A BITCH)
	data_labels = data_labels.replace(/'/g, '"')
	data_values = data_values.replace(/'/g, '"')

	var chart = new Chart(ctx,
	{
	  type: 'line',
	  data:
		{
	    labels: JSON.parse(data_labels),
	    datasets:
			[{
	        data: JSON.parse(data_values),
	        label: data_descr,
	        borderColor: "#009419",
	        fill: true
	     }]
	  },
	  options:
		{
	    title:
			{
	      display: true,
	      text: data_title
	    },
			legend:
			{
				display: false
			},
			scales:
			{
				yAxes:
				[{
					scaleLabel:
					{
		        display: true,
		        labelString: data_yaxis
					}
				}],
				xAxes:
				[{
					scaleLabel:
					{
		        display: true,
		        labelString: data_xaxis
					}
				}]
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

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
function barGraph(identifier, suffix)
{
	var target_id = identifier + suffix;
	var target_element = document.getElementById(target_id)
	// quit if null i.e. no data i.e. sample failed
	if (target_element === null){
		return;
	}
	var ctx = target_element.getContext('2d');
	var gradient = ctx.createLinearGradient(0,0,0,600);
	gradient.addColorStop(0, '#009419');
	gradient.addColorStop(1, '#002f08');

	// get data from JSON attribute
	var data_title = target_element.getAttribute('data-title')
	var data_label = target_element.getAttribute('data-descr')
	var data_labels = target_element.getAttribute('data-labels')
	var data_values = target_element.getAttribute('data-values')
	var data_xaxis = target_element.getAttribute('data-xaxis')
	var data_yaxis = target_element.getAttribute('data-yaxis')
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
						label: data_label,
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
					},
					scaleLabel:
					{
						display: true,
						labelString: data_yaxis
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
					},
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

// Test function to render dual barGraphs
function dualGraph(identifier, suffix)
{
	var target_id = identifier + suffix;
	var target_element = document.getElementById(target_id)
	// quit if null i.e. no data i.e. sample failed
	if (target_element === null){
		return;
	}
	var ctx = target_element.getContext('2d');
	var pri_gradient = ctx.createLinearGradient(0,0,0,600);
	pri_gradient.addColorStop(0, '#009419');
	pri_gradient.addColorStop(1, '#002f08');
	var sec_gradient = ctx.createLinearGradient(0,0,0,600);
	sec_gradient.addColorStop(0, '#94007b');
	sec_gradient.addColorStop(1, '#610051');

	// get data from JSON attribute
	var data_title = target_element.getAttribute('data-title')
	var pri_data_label = target_element.getAttribute('data-pri-descr')
	var pri_data_values = target_element.getAttribute('data-pri-values')
	var sec_data_label = target_element.getAttribute('data-sec-descr')
	var sec_data_values = target_element.getAttribute('data-sec-values')
	var data_labels = target_element.getAttribute('data-labels')
	var data_xaxis = target_element.getAttribute('data-xaxis')
	var data_yaxis = target_element.getAttribute('data-yaxis')
	// replace single quotes from python strings to double quotes (JSON A BITCH)
	data_labels = data_labels.replace(/'/g, '"')
	pri_data_values = pri_data_values.replace(/'/g, '"')
	sec_data_values = sec_data_values.replace(/'/g, '"')

	var chart = new Chart(ctx,
	{
		type: 'bar',
		data:
		{
				labels: JSON.parse(data_labels),
				datasets:
				[{
						label: pri_data_label,
						data: JSON.parse(pri_data_values),
						backgroundColor: pri_gradient,
						hoverBackgroundColor: pri_gradient,
						hoverBorderWidth: 2,
						hoverBorderColor: '#002f08'
				},
				{
					label: sec_data_label,
					data: JSON.parse(sec_data_values),
					backgroundColor: sec_gradient,
					hoverBackgroundColor: sec_gradient,
					hoverBorderWidth: 2,
					hoverBorderColor: '#610051'
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
			legend:
			{
				display: true
			},
			// Container for zoom options
      zoom:
			{
	      enabled: true,
	      mode: 'xy',
      },
			pan:
			{
				enabled: true,
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
					},
					scaleLabel:
					{
						display: true,
						labelString: data_yaxis
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
					},
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

// Function to render Box/Whisker/Line overlay graph
function boxGraph(identifier, suffix)
{
	var target_id = identifier + suffix;
	var target_element = document.getElementById(target_id)
	var ctx = target_element.getContext('2d');

	// data from attributes will go here
	var data_title = target_element.getAttribute("data-title")
	var data_labels = target_element.getAttribute("data-labels")
	var data_values = target_element.getAttribute("data-values")
	var data_meanval = target_element.getAttribute("data-meanval")
	var data_descr = target_element.getAttribute("data-descr")
	var data_xaxis = target_element.getAttribute("data-xaxis")
	var data_yaxis = target_element.getAttribute("data-yaxis")
	// replace single quotes from python strings to double quotes (JSON A BITCH)
	data_labels = data_labels.replace(/'/g, '"')
	data_values = data_values.replace(/'/g, '"')
	data_meanval = data_meanval.replace(/'/g, '"')

	var chart = new Chart(ctx,
		{
	    type: 'boxplot',
	    data:
			{
				labels: JSON.parse(data_labels),
				datasets: [{
					type: "boxplot",
					label: data_descr,
					backgroundColor: '#dcdcdc',
					borderColor: '#009419',
					borderWidth: 2,
					outlierColor: '#999999',
					padding: 10,
					itemRadius: 0,
					outlierColor: '#999999',
					data: JSON.parse(data_values)
				},
				{
					type: "line",
					label: "Mean PHRED score",
					borderColor: '#94007b',
					borderWidth: 2,
					fill: false,
					data: JSON.parse(data_meanval)
				}]
			},
	    options:
			{
	      responsive: true,
	      legend:
				{
					display: true
	      },
	      title:
				{
	        display: true,
	        text: data_title
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

// Function to render sequence alignment view
function sequenceView(identifier, allele_stage)
{
	var target_id = identifier + allele_stage;
	var target_element = document.getElementById(target_id);
	var input_string = target_element.getAttribute('data-sequences');

	// parsed array of the sequences
	var seqs =  msa.io.fasta.parse(input_string);
	var m = msa(
	{
     el: target_element,
     seqs: seqs,
		 bootstrapMenu: false,
		 colorscheme:
		 {
			 scheme: "nucleotide"
		 },
		 vis:
		 {
			 labels: false,
			 labelId: false,
			 overviewbox: false,
		 },
		 zoomer:
		 {
			 residueFont: "Product Sans",
			 alignmentHeight: 500,
		 },
		 conf:
		 {
			 debug: true
		 }

	});
	m.render();
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
		else {
			allParents[i].setAttribute('style', 'color: #009419;');
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
	$('html,body').animate({scrollTop: 0},'fast');

	// render Quality Control Graphs
	boxGraph(identifier, '_FQC_PBSQ');
	lineGraph(identifier, '_FQC_PBNC');
	lineGraph(identifier, '_FQC_SEQLENDIST');
	// render sequence alignment
	sequenceView(identifier, '_primary_seqmap');
	sequenceView(identifier, '_secondary_seqmap');
	// render Genotype Graphs
	barGraph(identifier, '_CCGDIST');
	dualGraph(identifier, '_CAGDIST');

	// Update header to be "SAMPLE | LANDING_STRING"
	rightCandidate = "\
	<p class=\"alignleft\">" + identifier + " | ScaleHD results" + "</p>\n\
	<p class=\"alignright\"><a class=\"contact\" href=\"mailto:alastair.maxwell@glasgow.ac.uk\"><strong>@</strong></a>  |  <a class=\"help\" href=\"#\" data-title=\"Help\"><strong>?</strong></a></p>\n\
	<div style=\"clear: both;\"></div>\n\
	"

	changeHeader(rightCandidate)

});

$('.sequence_sample_link_table').click(function(event){
	var identifier = $(event.target).data('sequenceid');
	var title = $(event.target).data('title');
	$('#help').hide();
	$('#welcome').hide();
	event.preventDefault();
	clearSampleColour(identifier)
	showSequence(identifier);
	$('html,body').animate({scrollTop: 0},'fast');

	//scroll the sidelist as well
	target_listIndex = $('#' + identifier + '_sampleLink').index() + 1;
	$('#sideNav').animate({scrollTop: $('#sideNav li:nth-child(' + target_listIndex + ')').position().top}, 'slow');

	// render Quality Control Graphs
	boxGraph(identifier, '_FQC_PBSQ');
	lineGraph(identifier, '_FQC_PBNC');
	lineGraph(identifier, '_FQC_SEQLENDIST');
	// render sequence alignment
	sequenceView(identifier, '_primary_seqmap');
	sequenceView(identifier, '_secondary_seqmap');
	// render Genotype Graphs
	barGraph(identifier, '_CCGDIST');
	dualGraph(identifier, '_CAGDIST');

	// Update header to be "SAMPLE | LANDING_STRING"
	rightCandidate = "\
	<p class=\"alignleft\">" + identifier + " | ScaleHD results" + "</p>\n\
	<p class=\"alignright\"><a class=\"contact\" href=\"mailto:alastair.maxwell@glasgow.ac.uk\"><strong>@</strong></a>  |  <a class=\"help\" href=\"#\" data-title=\"Help\"><strong>?</strong></a></p>\n\
	<div style=\"clear: both;\"></div>\n\
	"

	changeHeader(rightCandidate)

});

// When the #home div button is pressed
$('.home').click(function(event){
	clearSampleColour('None')
	scrollToTop();
	$('div.sequence').hide();
	$('#help').hide();
	$('#welcome').show();
	barGraph(instanceLabel, 'CAGSummaryChart');
	barGraph(instanceLabel, 'CCGSummaryChart');

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
	// render quality control graphs
	boxGraph(identifier, '_FQC_PBSQ');
	lineGraph(identifier, '_FQC_PBNC');
	lineGraph(identifier, '_FQC_SEQLENDIST');
	// render sequence alignment
	sequenceView(identifier, '_primary_seqmap');
	sequenceView(identifier, '_secondary_seqmap');
	// render genotype graphs
	barGraph(identifier, '_CCGDIST');
	dualGraph(identifier, '_CAGDIST');

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

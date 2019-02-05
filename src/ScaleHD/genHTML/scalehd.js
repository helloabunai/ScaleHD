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

	console.log(target_element.attr('data-test'))

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
function CAGSummaryChart(identifier){
	var target = identifier + "CAGSummaryChart";
	var ctx = document.getElementById(target).getContext('2d');
	var gradient = ctx.createLinearGradient(0,0,0,600);
	gradient.addColorStop(0, '#009419');
	gradient.addColorStop(1, '#001504');
	var myChart = new Chart(ctx,{
	    type: 'bar',
	    data:
			{
	        labels: ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '200'],
	        datasets:
					[{
	            label: '# of alleles present',
	            data: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	            backgroundColor: gradient,
							hoverBackgroundColor: gradient,
							hoverBorderWidth: 2,
							hoverBorderColor: '#001504'
	        }]
	    },
	    options:
			{
				// Prevents graphs resizing
				responsive: false,
				// Title colour/text
				title:
				{
            display: true,
						fontColor: '#000000',
						fontSize: 15,
						fontFamily: 'Product Sans',
            text: "CAG allele distribution for webtest"
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


// CCG Summary
function CCGSummaryChart(identifier){
	var target = identifier + "CCGSummaryChart";
	var ctx = document.getElementById(target).getContext('2d');
	var gradient = ctx.createLinearGradient(0,0,0,600);
	gradient.addColorStop(0, '#009419');
	gradient.addColorStop(1, '#001504');
	var myChart = new Chart(ctx,{
	    type: 'bar',
	    data:
			{
	        labels: ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'],
	        datasets:
					[{
	            label: '# of alleles present',
	            data: [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	            backgroundColor: gradient,
							hoverBackgroundColor: gradient,
							hoverBorderWidth: 2,
							hoverBorderColor: '#001504'
	        }]
	    },
	    options:
			{
				// Prevents graphs resizing
				responsive: false,
				// Title colour/text
				title:
				{
            display: true,
						fontColor: '#000000',
						fontSize: 15,
						fontFamily: 'Product Sans',
            text: "CCG allele distribution for webtest"
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


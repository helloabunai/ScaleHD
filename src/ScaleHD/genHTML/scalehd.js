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
	summaryCAGChart(instanceLabel);
	summaryCCGChart(instanceLabel);
});

// Function for clicking on a sequence sample
function showSequence(identifier){
	$('div.sequence').hide();
	$('#welcome').hide();
	var seqHead = $('div.sequence[data-sequenceid='+ identifier +']')
	seqHead.show();
	//$('html,body').animate({scrollTop: seqHead.offset().top - 50}, 'fast')
}

// Function to change HTML header text if sample/substage link pressed
function changeHeader(replaceRight){
	$('#rightHeader').html(replaceRight)
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

////////////////
// Chart Test //
////////////////

// CAG Summary
function summaryCAGChart(identifier)
{
	var target = identifier + 'CAGSummaryChart';
	var ctx = document.getElementById(target).getContext('2d');
	var myChart = new Chart(ctx, {
	    type: 'bar',
	    data: {
	        labels: ["Red", "Blue", "Yellow", "Green", "Purple", "Orange"],
	        datasets: [{
	            label: '# of alleles',
	            data: [12, 19, 3, 5, 2, 3],
	            backgroundColor: [
	                'rgba(255, 99, 132, 0.2)',
	                'rgba(54, 162, 235, 0.2)',
	                'rgba(255, 206, 86, 0.2)',
	                'rgba(75, 192, 192, 0.2)',
	                'rgba(153, 102, 255, 0.2)',
	                'rgba(255, 159, 64, 0.2)'
	            ],
	            borderColor: [
	                'rgba(255,99,132,1)',
	                'rgba(54, 162, 235, 1)',
	                'rgba(255, 206, 86, 1)',
	                'rgba(75, 192, 192, 1)',
	                'rgba(153, 102, 255, 1)',
	                'rgba(255, 159, 64, 1)'
	            ],
	            borderWidth: 1
	        }]
	    },
	    options: {
				responsive: false,
				title: {
            display: true,
            text: 'CAG test'
        },
	      scales: {
	   			yAxes: [{
	        	ticks: {
	          	beginAtZero:true
	            }
	           }]
	        }
	    }
	});
}

// CCG Summary
function summaryCCGChart(identifier)
{
	var target = identifier + 'CCGSummaryChart';
	var ctx = document.getElementById(target).getContext('2d');
	var myChart = new Chart(ctx, {
	    type: 'bar',
	    data: {
	        labels: ["Red", "Blue", "Yellow", "Green", "Purple", "Orange", "test1", "test2", "test3", "test4", "test5"],
	        datasets: [{
	            label: '# of alleles',
	            data: [12, 19, 3, 5, 2, 3, 20, 14, 9, 8, 12],
	            backgroundColor: [
	                'rgba(255, 99, 132, 0.2)',
	                'rgba(54, 162, 235, 0.2)',
	                'rgba(255, 206, 86, 0.2)',
	                'rgba(75, 192, 192, 0.2)',
	                'rgba(153, 102, 255, 0.2)',
	                'rgba(255, 159, 64, 0.2)',
									'#666666',
									'#666666',
									'#666666',
									'#666666',
									'#666666'
	            ],
	            borderColor: [
	                'rgba(255,99,132,1)',
	                'rgba(54, 162, 235, 1)',
	                'rgba(255, 206, 86, 1)',
	                'rgba(75, 192, 192, 1)',
	                'rgba(153, 102, 255, 1)',
	                'rgba(255, 159, 64, 1)',
									'#666666',
									'#666666',
									'#666666',
									'#666666',
									'#666666'
	            ],
	            borderWidth: 1
	        }]
	    },
			options: {
				responsive: false,
				title: {
            display: true,
            text: 'CCG test'
        },
	      scales: {
	   			yAxes: [{
	        	ticks: {
	          	beginAtZero:true
	            }
	           }]
	        }
	    }
	});
}

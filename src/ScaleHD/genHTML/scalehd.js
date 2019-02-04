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

	TEST_BOX(identifier);

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


// test box Plot
function randomValues(count, min, max) {
  const delta = max - min;
  return Array.from({length: count}).map(() => Math.random() * delta + min);
}

const boxplotData = {
  // define label tree
  labels: ['January', 'February', 'March', 'April', 'May', 'June', 'July'],
  datasets: [{
    label: 'Dataset 1',
    backgroundColor: 'rgba(255,0,0,0.5)',
    borderColor: 'red',
    borderWidth: 1,
    outlierColor: '#999999',
    padding: 10,
    itemRadius: 0,
    outlierColor: '#999999',
    data: [
      randomValues(100, 0, 100),
      randomValues(100, 0, 20),
      randomValues(100, 20, 70),
      randomValues(100, 60, 100),
      randomValues(40, 50, 100),
      randomValues(100, 60, 120),
      randomValues(100, 80, 100)
    ]
  }, {
    label: 'Dataset 2',
    backgroundColor:  'rgba(0,0,255,0.5)',
    borderColor: 'blue',
    borderWidth: 1,
    outlierColor: '#999999',
    padding: 10,
    itemRadius: 0,
    outlierColor: '#999999',
    data: [
      randomValues(100, 60, 100),
      randomValues(100, 0, 100),
      randomValues(100, 0, 20),
      randomValues(100, 20, 70),
      randomValues(40, 60, 120),
      randomValues(100, 20, 100),
      randomValues(100, 80, 100)
    ]
  }]
};

function TEST_BOX(identifier){
	var target = "FQC_PBSQ_" + identifier;
	console.log(target)
	var ctx = document.getElementById(target).getContext('2d');
	var myChart = new Chart(ctx, {
    type: 'boxplot',
    data: boxplotData,
    options: {
      responsive: true,
      legend: {
        position: 'top',
      },
      title: {
        display: true,
        text: 'Chart.js Box Plot Chart'
      }
    }
  });

};

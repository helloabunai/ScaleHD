function {TRIPLET_TYPE}(identifier)
{
	var target = identifier + "{TRIPLET_TYPE}";
	var ctx = document.getElementById(target).getContext('2d');
	var myChart = new Chart(ctx, {
	    type: 'bar',
	    data: {
	        labels: {LABELS},
	        datasets: [{
	            label: '# of alleles present',
	            data: {TRIPLET_DISTRIBUTION},
	            fillColor: {BACKGROUND_COLOUR},
	            strokeColor:{BORDER_COLOUR},
	            borderWidth: 1
	        }]
	    },
	    options: {
				responsive: false,
				title: {
            display: true,
            text: {CHART_TITLE}
        },
				legend:{
					display: false
				},
	      scales: {
	   			yAxes: [{
	        	ticks: {
	          	beginAtZero:true,
							precision:0
	            }
						}],
					xAxes: [{
						ticks: {
							max:200,
							min:0,
							stepSize: 20
							}
						}]
	        }
	    }
	});
}

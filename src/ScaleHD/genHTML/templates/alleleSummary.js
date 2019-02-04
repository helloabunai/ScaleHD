function {TRIPLET_TYPE}(identifier){
	var target = identifier + "{TRIPLET_TYPE}";
	var ctx = document.getElementById(target).getContext('2d');
	var gradient = ctx.createLinearGradient(0,0,0,600);
	gradient.addColorStop(0, '#009419');
	gradient.addColorStop(1, '#001504');
	var myChart = new Chart(ctx,{
	    type: 'bar',
	    data:
			{
	        labels: {LABELS},
	        datasets:
					[{
	            label: '# of alleles present',
	            data: {TRIPLET_DISTRIBUTION},
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
            text: {CHART_TITLE}
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

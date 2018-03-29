function click_dummy_button(){
	// this click the hidden button in "dum_box" in order to start/stop the counter in the status_div
	document.getElementsByClassName("dum_button")[0].children[0].click();
}

/*
document.addEventListener("DOMContentLoaded", function(event){
	setTimeout( function(){
		console.log(document)

		//document.getElementsByClassName('main_grid')[0].style["width"] = "1215px"; // this is to allow the addition of a border to 'side_box'
		document.getElementsByClassName('custom_button')[0].parentElement.style["padding"] = "0px";

	}, 2000)
})
*/
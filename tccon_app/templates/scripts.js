function click_dummy_button(){
	// this click the button in the hidden "dum_box" in order to start/stop the counter in the status_div
	document.getElementsByClassName("dum_button")[0].children[0].click();
}

document.addEventListener("DOMContentLoaded", function(event){
	setTimeout( function(){
		console.log(document)

		//seems like setting those in styles.css does not work.
		document.getElementsByClassName('dum_box')[0].style["display"] = "none"; 
		document.getElementsByClassName('main_grid')[0].style["width"] = "1215px";
		document.getElementsByClassName('custom_button')[0].parentElement.style["padding"] = "0px";

	}, 2000)
})
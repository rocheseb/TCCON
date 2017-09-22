function click_dummy_button(){
	var dum_button = document.getElementsByClassName("dum_button")[0].children[0];
	dum_button.click();
}

document.addEventListener("DOMContentLoaded", function(event){
	setTimeout( function(){
		console.log(document)
		var dum_box = document.getElementsByClassName("dum_button")[0].parentElement;
		dum_box.style.display="none";

		var toolbar_list = document.getElementsByClassName('bk-toolbar-box');
		if (toolbar_list.length > 1){
			for(var i=0;i < toolbar_list.length;i++){
				if (toolbar_list[i].textContent.includes('Hover')){
					toolbar_list[i].style.display = 'none'
				}
			}
		}
	}, 2000)
})
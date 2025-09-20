

 /* src/05b_e_mediarecorder.js - start*/ 
 const registra = {};


 registra.template = `<div id="videoContainer">
                       <div id="videoPlayerBox">
                       <div class="closeVideoBox">x</div>
                       </div>
                     </div>`;
 
 registra.start = function(){
   console.log('>>>>>  registra.start ');
   let recording = false;
   let chunks = [];
   
   const framerate_const = 30;
   
   function record() {
       chunks.length = 0;
       let stream = document.querySelector('canvas').captureStream(framerate_const);
       registra.recorder = new MediaRecorder(stream );
       registra.recorder.ondataavailable = e => {
           if (e.data.size) {
               chunks.push(e.data);
           }
       };
       registra.recorder.onstop = exportVideo;
   }
   
   function exportVideo(e) {
     let fileName = prompt("Please enter your video clip name");
     if(fileName && fileName !== null && fileName !== "null") {
       var blob = new Blob(chunks, { 'type': 'video/webm' });
       // Draw video to screen
       $('body').append(registra.template);
       var videoElement = document.createElement('video');
       videoElement.setAttribute("id", Date.now());
       videoElement.controls = true;
       $('#videoPlayerBox').append(videoElement);
       videoElement.src = window.URL.createObjectURL(blob);
 
       // Download the video 
       var url = URL.createObjectURL(blob);
       var a = document.createElement('a');
       a.setAttribute("id", "downloadVideo");
       a.innerHTML = "Download video";
       a.style.display = 'none';
       $('#videoPlayerBox').append(a);
       a.href = url;
       a.download = fileName+'.webm';
       a.click();
       window.URL.revokeObjectURL(url);
 
       }else{
         console.log('>>>>> fileName undefined or null')
       }
   }
   record();
 }
 
 
 /*
 function keyPressed() {
     // toggle recording true or false
     recording = !recording
     console.log(recording);
     // 82 is keyCode for r 
     // if recording now true, start recording 
     if (keyCode === 82 && recording ) {
       console.log("recording started!");
       registra.recorder.start();
     } 
     // if we are recording, stop recording 
     if (keyCode === 82 && !recording) {  
       console.log("recording stopped!");
       registra.recorder.stop();
     }
   }
   */
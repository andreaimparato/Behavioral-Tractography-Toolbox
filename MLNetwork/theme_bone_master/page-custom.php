<?php
/*
 Template Name: 3dVisualizer HP
* Template Post Type: page
 *
 * This is your custom page template. You can create as many of these as you need.
 * Simply name is "page-whatever.php" and in add the "Template Name" title at the
 * top, the same way it is here.
 *
 * When you create your page, you can just select the template and viola, you have
 * a custom page template to call your very own. Your mother would be so proud.
 *
 * For more info: http://codex.wordpress.org/Page_Templates
*/
?>


<!doctype html>

<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">

  <title>MLNETWORK-DIPLAB-CH</title>
  <meta name="description" content="MLNETWORK-DIPLAB-CH">
  <meta name="author" content="">

  <meta property="og:title" content="MLNETWORK-DIPLAB-CH">
  <meta property="og:type" content="website">
  <meta property="og:url" content="">
  <meta property="og:description" content="">
  <meta property="og:image" content="image.png">

  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
  <link rel="icon" href="<?php echo get_template_directory_uri(); ?>/favicon.ico">
  <link rel="icon" href="/favicon.svg" type="image/svg+xml">
  <link rel="apple-touch-icon" href="/apple-touch-icon.png">

  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Poppins:ital,wght@0,100;0,200;0,300;0,400;0,500;0,600;0,700;0,800;0,900;1,100;1,200;1,300;1,400;1,500;1,600;1,700;1,800;1,900&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="<?php echo get_template_directory_uri(); ?>/css/style.css?v=3">
  <link rel="stylesheet" href="<?php echo get_template_directory_uri(); ?>/css/style1.css?v=3">
  <!-- PLEASE NO CHANGES BELOW THIS LINE (UNTIL I SAY SO) -->
  <script language="javascript" type="text/javascript" src="https://cdn.jsdelivr.net/npm/lodash@4.17.21/lodash.min.js"></script>
  <script language="javascript" type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/p5.js/1.6.0/p5.min.js"></script>
  <script language="javascript" src="<?php echo get_template_directory_uri(); ?>/_library/p5.easycam.js"></script>
  <script language="javascript" type="text/javascript" src="<?php echo get_template_directory_uri(); ?>/dist/app.min.js?v=4"></script>

</head>

<body>

<?php if (have_posts()) : while (have_posts()) : the_post(); ?>



    
    <!-- MENU -->
    <div class="menu">
        <div class="menu-collapse"><button></button></div>
        <div class="logo">Multi-Layer Network Analyzer</span></div>
        <div class="comands">
            <div class="title"><span>Commands</span></div>
            <div class="comand-items">
                <div class="comand-item">
                    <div class="label"> Distance between layers</div>
                    <div class="input-minus-plus layers">
                        <button class="meno">-</button>
                        <button class="piu">+</button>
                    </div>
                </div>
                <div class="comand-item">
                    <div class="label">Node size proportion</div>
                    <div class="input-minus-plus nodi">
                        <button class="meno">-</button>
                        <button class="piu">+</button>
                    </div>
                </div>


                <div class="comand-item">
                    <div class="label">Edge width proportion</div>
                    <div class="input-minus-plus edge">
                        <button class="meno">-</button>
                        <button class="piu">+</button>
                    </div>
                </div>


                <div class="comand-item comand-item-button layer">
                    <div class="label button">Show/Hide Layers</div>
                </div>
                <div class="menu-layer listalayer">
                </div>


                <div class="comand-item comand-item-button nodi">
                    <div class="label button">Select Nodes</div>
                </div>
                <div class="menu-layer listanodi">
                    <div class="nodi_wrapper">
                    </div>
                </div>

            </div>
        </div>
    </div>
    <!-- END - MENU -->

    <!-- DRAWER NODI -->
    <div class="drawer drawer_nodi">
        <div class="wrapper">
            <div class="close"><a href="#">close</a></div>
            <div class="drawer-title">Select Nodes</div>
            <div class="dynamic-content">
                <div class="load-message">Load JSON</div>
            </div>
        </div>
    </div>
    <!-- END - DRAWER NODI -->

    <!-- DRAWER LAYERS -->
    <div class="drawer drawer_layers">
        <div class="wrapper">
            <div class="close"><a href="#">close</a></div>
            <div class="drawer-title">Select Layers</div>
            <div class="dynamic-content">
                <div class="load-message">Load JSON</div>
            </div>
        </div>
    </div>
    <!-- END - DRAWER LAYERS -->



    <div class="comandi-bottom">
        <div class="comandi-wrapper">
            <div class="comandi-collapse"><button></button></div>
            <div class="group">
                <div class="comando-bottom-item-title">Layout</div>
                <div class="item">
                    <div class="item-title">Background</div>
                    <div class="item-input">
                        <input type="color" value="#129b99" id="canvasColor">
                    </div>
                </div>
            </div>
            <div class="group">
                <div class="comando-bottom-item-title">Movement</div>
                <div class="item">
                    <div class="item-title">Graph</div>
                    <div class="item-input input-arrows">
                        <button class="arrow arrow-up" onclick="easycam.rotateX(0.1)"></button>
                        <button class="arrow arrow-down" onclick="easycam.rotateX(-0.1)"></button>
                        <button class="arrow arrow-left" onclick="easycam.rotateY(0.1)"></button>
                        <button class="arrow arrow-right" onclick="easycam.rotateY(-0.1)"></button>
                    </div>
                </div>
                <div class="item">
                    <div class="item-title">Label</div>
                    <div class="item-input input-arrows">
                        <button class="arrow arrow-up" onclick="textRotationY(+.1)"></button>
                        <button class="arrow arrow-down" onclick="textRotationY(-.1)"></button>
                        <button class="arrow arrow-left" onclick="textRotationX(+.1)"></button>
                        <button class="arrow arrow-right" onclick="textRotationX(-.1)"></button>
                    </div>
                </div>

                <div class="comando-bottom-item-title labelSize">
                    <div class="labelDecSize" onclick="labelSetSize('-')">-</div>
                    <div class="labelIncSize" onclick="labelSetSize('+')">+</div>
                </div>
            </div>
        </div>  
    </div>


    <div class="textPopup">
        <div class="labelProperty">
            <div class="fontsize"></div>
            <div class="color"></div>
            <div class="position position_contaniner" xxx-nodo-id="">
                <div class="position_item">
                    <div class="position_item_title">Display Name</div>
                    <div class="item_name"><input type="text" name="name" value=""></div>
                </div>
                <div class="position_item">
                    <div class="position_item_title">Text Style</div>
                    <div class="position_item_controls">
                        <span class="control_button regular style_button" xxx-data-style="regular">Regular</span>
                        <span class="control_button bold style_button" xxx-data-style="bold">Bold</span>
                    </div>
                </div>

                <div class="position_item">
                    <div class="position_item_title">Label size</div>
                    <div class="position_item_controls labelSize">
                        <span class="control_button minus style_button" xxx-node-index="" onclick="labelSetSizeSingular(this.getAttribute('xxx-data-style'), this.getAttribute('xxx-node-index'))" xxx-data-style="-">-</span>
                        <span class="control_button plus style_button" xxx-node-index="" onclick="labelSetSizeSingular(this.getAttribute('xxx-data-style'), this.getAttribute('xxx-node-index') )" xxx-data-style="+">+</span>
                    </div>
                </div>
                <div class="position_item">
                    <div class="position_item_title">Position</div>
                    <div class="position_item_controls">
                        <span class="control_button top position_button" xxx-data-position="top">top</span>
                        <span class="control_button bottom position_button" xxx-data-position="bottom">bottom</span>
                    </div>
                </div>
                <div class="position_item">
                    <div class="position_item_title">Translate</div>
                    <div class="position_item_controls">
                        <span class="left_num" >
                            <input id="labelAdjustY" name="nodo_left">
                            <label>X AXIS</label>
                        </span>
                        <span class="right_num">
                            <input id="labelAdjustX" name="nodo_right">
                            <label>Y AXIS</label>
                        </span>
                    </div>
                </div>
                <div class="position_item">
                    <div class="position_item_title">Color</div>
                    <div class="position_item_controls">
                        <input type="color" value="#000000" id="nodeColor" onchange="nodeColorChange(this)">
                    </div>
                </div>

                <div class="button">
                    <div class="buttonCancel" onclick="closeLabelPopup()">Cancel</div>
                    <div class="buttonSave" onclick="saveLabelPopup()">Save</div>
                </div>
            </div>
        </div>
    </div>

    <!-- BODY -->
    <div class="body">
        <div class="header">
            <div class="intro"></div>
            <div class="interactions">
                <div class="load-file">
                    <!--<div class="file-loaded">File caricato<span>nome-file.json</span></div>-->
                    <form id="jsonFile" name="jsonFile" enctype="multipart/form-data" method="post">
                        <fieldset>
                           <input type='file' id='fileinput'>
                        </fieldset>
                      </form>
                    <!--<button>Scegli file</button>-->
                </div>
                <div class="save-image">
                    <button onclick="savec()"><span>Save image</span></button>
                </div>
                <div class="download-json">
                    <button onclick="exportFile()"><span>Save json</span></button>
                </div>
                <div class="screen-recording" xxx-video-recorder="false">
                    <button class="start-recording" onclick="registra.recorder.start()"><span class="rec-led"></span><span>Start REC</span></button>
                    <button class="stop-recording" onclick="registra.recorder.stop()"><span class="rec-led"></span><span>Stop REC</span></button>
                </div>
            </div>
        </div>
    </div>

    <!-- END - BODY -->


	<?php endwhile; endif; ?>
</body>
</html>



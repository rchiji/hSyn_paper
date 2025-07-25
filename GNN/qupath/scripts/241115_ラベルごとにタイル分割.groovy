List annotations = getAnnotationObjects().findAll {
    it.getPathClass().toString() != "Region*"
    }

Set newTiles = []

annotations.each { anno -> 
    def parentROI = anno.getROI()
    List tileROIs = RoiTools.makeTiles(parentROI,50,50,true)  
    
    tileROIs.each { tileROI -> 
        if ( tileROI.getArea() >= 20 ) {
            newTiles.add( PathObjects.createTileObject( tileROI, anno.getPathClass() ) )
            } // <- if ( tileROI.getArea() >= 50 )

        } // <- tileROIs.each   
    } // <- annotations.each

addObjects(newTiles)

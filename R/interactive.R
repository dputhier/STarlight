#' Interactive Visual Feature Picker
#'
#' This function launches an interactive Shiny application for visual exploration
#' and selection of features within a given dataset. It allows users to interactively
#' select features based on various criteria and visualization tools such as tangents,
#' areas, nearest points, and polygons.
#'
#' @param object A dataframe or a matrix that contains the coordinates (x, y) and
#'               the feature classification. It should have at least three columns:
#'               x (x-coordinates), y (y-coordinates), and feature (feature labels).
#'               Additional columns are supported.
#'
#' @return Launches a Shiny application which does not return a value to the R environment.
#'         Selected features and configurations can be exported directly from the application interface.
#'
#' @examples
#' #' \dontrun{
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' visual_feature_picker(coord(xen))
#' data_set <- data.frame(x = runif(1000, 0, 1000), y = runif(1000, 0, 1000), feature = sample(c("Feature1", "Feature2"), 100, replace = TRUE))
#' visual_feature_picker(data_set)
#' }
#' @importFrom shiny reactiveValues shinyApp selectInput actionButton renderUI
#'              uiOutput selectizeInput renderPlot observeEvent downloadHandler
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#'              dashboardBody sidebarMenu menuItem box fluidRow
#' @importFrom ggplot2 ggplot aes geom_point theme_bw scale_color_manual theme
#'              element_text element_blank guide_legend geom_abline geom_segment
#' @importFrom colourpicker colourInput
#' @importFrom secr pointsInPolygon
#' @export
#'
visual_feature_picker <- function(object=NULL) {
  if(is.null(object))
    print_this_msg("Need an input object.")
  initialize_data <- function() {
    tmp <- object[, c("x", "y", "feature")]
    tmp$color <- "selected"
    tmp$to_be_checked <- 0

    ## Set up the reactive dataframe
    my_val <- shiny::reactiveValues()
    my_val$DT <- tmp
    my_val$gline <- NULL
    my_val$gpoint <- NULL
    my_val$ggseg <- NULL

    return(my_val)
  }

  draw_plot <- function(input, my_val) {
    print_this_msg("Representative feature is:", input$repr_feature)
    # Avoids some warnings
    if (all(c('selected', 'unselected') %in% my_val$DT$color)) {
      print_this_msg("The dataframe contains selected/unselected values.",
                     msg_type = "DEBUG")
      scale_prep <- ggplot2::scale_color_manual(values = c(
        'selected' = input$selected_col,
        'unselected' = input$unselected_col
      ))
    } else if ('selected' %in% my_val$DT$color &
               !'unselected' %in% my_val$DT$color) {
      print_this_msg("The dataframe contains selected values only.", msg_type = "DEBUG")
      scale_prep <-
        ggplot2::scale_color_manual(values = c('selected' = input$selected_col))
    } else if (!'selected' %in% my_val$DT$color &
               'unselected' %in% my_val$DT$color) {
      print_this_msg("The dataframe contains unselected values only.",
                     msg_type = "DEBUG")
      scale_prep <-
        ggplot2::scale_color_manual(values = c('unselected' = input$unselected_col))
    }

    print_this_msg("Preparing ggplot", msg_type = "DEBUG")
    ggprep <-
      ggplot2::ggplot(my_val$DT[my_val$DT$feature == input$repr_feature, ],
                      ggplot2::aes(x = x, y = y, color = color)) +
      ggplot2::geom_point(size = as.double(input$pts_size))  +
      ggplot2::theme_bw() +
      scale_prep +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_blank(),
        legend.position = "bottom"
      ) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size =
                                                                           4,
                                                                         shape =
                                                                           15))) +
      ggplot2::coord_fixed()


    if (input$tool == "Tangents") {
      print_this_msg("Inside drawplot -> Tangents.", msg_type = "DEBUG")
      if (is.null(my_val$gline)) {
        gline <- NULL
      } else{
        gline <- my_val$gline
      }

      if (is.null(my_val$gpoint)) {
        gpoint <- NULL
      } else{
        gpoint <- my_val$gpoint
      }


      gg <- ggprep +
        gline +
        gpoint

    } else if (input$tool == "Areas") {
      print_this_msg("Inside drawplot -> Areas.", msg_type = "DEBUG")

      gg <- ggprep

    } else if (input$tool == "Nearest Points") {
      print_this_msg("Inside drawplot -> Nearest Points.", msg_type = "DEBUG")
      gg <- ggprep
    } else if (input$tool == "Polygons") {
      print_this_msg("Inside drawplot -> Polygons.", msg_type = "DEBUG")
      gg <- ggprep + my_val$ggseg
    }

    print_this_msg("Returning ggplot", msg_type = "DEBUG")

    return(gg)

  }

  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "Select features"),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Menu",
                                 tabName = "main",
                                 icon = shiny::icon("th")),
        shiny::uiOutput("feature_select")
        ,
        shiny::selectInput(
          "tool",
          "Tool:",
          choices = c("Tangents", "Areas", "Nearest Points", "Polygons")
        )
      ),
      shiny::actionButton("swap_points",
                          "Swap selected/unselected",
                          class = "btn-primary"),
      shiny::actionButton("reset",
                          "Reset to default",
                          class = "btn-primary"),
      shiny::actionButton("unselect_all",
                          "Unselect all",
                          class = "btn-primary"),
      shinydashboard::sidebarMenu(
        colourpicker::colourInput("selected_col", "Color for selected points", "black"),
        colourpicker::colourInput("unselected_col", "Color for unselected points", "lightgray"),
        shiny::textInput(
          "pts_size",
          "Point size:",
          value = "0.5",
          placeholder = "0.5"
        ),
        shiny::checkboxInput(
          "save_selected",
          "Save selected (default) or unselected points",
          TRUE
        ),
        shiny::downloadButton('download_data', 'Download', class = "dload"),
        htmltools::tags$head(
          htmltools::tags$style(
            ".dload{color: black !important; background-color: white !important}"
          )
        )
      )
    ),
    shinydashboard::dashboardBody(shinydashboard::tabItems(
      shinydashboard::tabItem(tabName = "main",
                              shiny::fluidRow(
                                shinydashboard::box(
                                  title = "Tool Description and parameters",
                                  status = "primary",
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  width = 12,
                                  shiny::uiOutput("outils_ui")
                                )
                              ),
                              shiny::fluidRow(
                                shinydashboard::box(
                                  title = "Main plot",
                                  status = "warning",
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  width = 12,
                                  height = "1000px",
                                  shiny::plotOutput("main_plot",
                                             click = "click",
                                             brush = "brush_select")
                                )
                              ))
    ))
  )

  server <- function(input, output, session) {
    my_val <- initialize_data()

    output$feature_select <- renderUI({
      shiny::selectizeInput(
        "repr_feature",
        "Representative feature",
        sort(unique(my_val$DT$feature)),
        selected = sort(unique(my_val$DT$feature))[1],
        multiple = FALSE,
        options = NULL
      )
    })

    ## Observe event (rem_points_tangents)
    shiny::observeEvent(input$repr_feature, {
      output$feature_select <- renderUI({
        shiny::selectizeInput(
          "repr_feature",
          "Representative feature",
          sort(unique(my_val$DT$feature)),
          selected = input$repr_feature,
          multiple = FALSE,
          options = NULL
        )
      })
    })

    output$outils_ui <- renderUI({
      switch(
        input$tool,
        "Tangents" = {
          htmltools::tagList(
            htmltools::h4("Select two points to create a straight line."),
            shiny::checkboxInput(
              "un_select_below",
              "Unselect point below (default) or above the line.",
              TRUE
            ),
            shiny::actionButton("rem_points_tangents", "Remove the points", class = "btn-primary")
          )
        },
        "Areas" = {
          htmltools::tagList(
            htmltools::h4("Click and drag to create areas."),
            shiny::checkboxInput(
              "area_unselect",
              "Areas are used to unselect (default) or to select points.",
              TRUE
            ),
          )
        },
        "Nearest Points" = {
          htmltools::tagList(
            htmltools::h4("Click to select the Nearest Points."),
            shiny::fluidRow(column(
              4,
              shiny::checkboxInput(
                "circular_unselect",
                "Circular area are used to unselect (default) or to select points",
                TRUE
              )
            )),
            shiny::fluidRow(column(
              2, shiny::textInput("radius", "Circle radius:", value = "500")
            ),
            column(
              2, shiny::textInput("threshold", "Max distance:", value = "5")
            ))
          )
        },
        "Polygons" = {
          htmltools::tagList(
            htmltools::h4("Select point one by one and click on 'Close the path'."),
            shiny::actionButton("close_the_path", "Close the path", class = "btn-primary")
          )
        }
      )
    })


    output$main_plot <- shiny::renderPlot({
      gp <- draw_plot(input, my_val)
      gp

    },
    height = 800,
    width = 800)

    ## Observe event (click)
    shiny::observeEvent(input$click, {
      if (input$tool == "Tangents") {
        print_this_msg("Inside input$click -> Tangents",
                       msg_type = "DEBUG")
        if (nrow(my_val$DT[my_val$DT$feature == "POINT_COORDS", ]) < 2) {
          add_row <- data.frame(
            x = input$click$x,
            y = input$click$y,
            feature = "POINT_COORDS",
            color = NA,
            to_be_checked = 0
          )

          my_val$DT <- rbind(my_val$DT, add_row)

          l_pts <- my_val$DT[my_val$DT$feature == "POINT_COORDS",]

          gpoint <- ggplot2::geom_point(data = l_pts,
                                        size = 4,
                                        color = "red")
          my_val$gpoint <- gpoint

          if (nrow(l_pts) == 2) {
            # calculating the slope
            a <-
              (l_pts[1, 2] - l_pts[2, 2]) / (l_pts[1, 1] - l_pts[2, 1])
            # calculating intersect
            b <- l_pts[1, 2] - l_pts[1, 1] * a
            df <- data.frame(cbind(a, b))
            gline <- ggplot2::geom_abline(data = df,
                                          ggplot2::aes(slope = a,
                                                       intercept = b))
            if (input$un_select_below) {
              my_val$DT[a * my_val$DT$x + b > my_val$DT$y, "color"] <-
                "unselected"
            } else{
              my_val$DT[a * my_val$DT$x + b < my_val$DT$y, "color"] <-
                "unselected"
            }

            my_val$gline <- gline
            my_val$DT <-
              my_val$DT[my_val$DT$feature != "POINT_COORDS", ]
          }

        }
      } else if (input$tool == "Nearest Points") {
        print_this_msg("Inside input$click -> Nearest Points",
                       msg_type = "DEBUG")
        shiny::nearPoints(
          my_val$DT,
          input$click,
          xvar = "x",
          yvar = "y",
          threshold = input$threshold
        )
        xc <- input$click$x
        yc <- input$click$y
        radius <- as.double(input$radius)
        inside <-
          sqrt((my_val$DT$x - xc) ^ 2 + (my_val$DT$y - yc) ^ 2) < radius

        print_this_msg("Number of points inside circle:",
                       sum(inside),
                       msg_type = "DEBUG")

        if (input$circular_unselect) {
          my_val$DT$color[inside] <- "unselected"
        } else{
          my_val$DT$color[inside] <- "selected"
        }


      } else if (input$tool == "Polygons") {
        print_this_msg("Inside input$click -> Polygons",
                       msg_type = "DEBUG")

        add_row <- data.frame(
          x = input$click$x,
          y = input$click$y,
          feature = "POLYGON_COORD",
          color = "Black",
          to_be_checked = 0
        )

        my_val$DT <- rbind(my_val$DT, add_row)

        seg_pts <-
          my_val$DT[my_val$DT$feature == "POLYGON_COORD",]

        if (nrow(seg_pts) >= 2) {
          seg_pts$xend <- NA
          seg_pts$yend <- NA

          for (i in 1:(nrow(seg_pts) - 1)) {
            seg_pts[i,]$xend <- seg_pts[i + 1,]$x
            seg_pts[i,]$yend <- seg_pts[i + 1,]$y
          }

          my_val$polygon <- seg_pts

          seg_pts <- na.omit(seg_pts)

          ggseg <- ggplot2::geom_segment(
            data = seg_pts,
            size = 1,
            color = "blue",
            aes(
              x = x,
              y = y,
              xend = xend,
              yend = yend
            )
          )

          my_val$ggseg <- ggseg
        }

      }

      output$main_plot <- shiny::renderPlot({
        draw_plot(input, my_val)
      },
      height = 800,
      width = 800)
    })

    ## Observe event (brush_select)
    shiny::observeEvent(input$brush_select, {
      print_this_msg("Inside input$brush_select event", msg_type = "DEBUG")
      print_this_msg("Tool:", input$tool, msg_type = "DEBUG")

      if (input$tool == "Areas") {
        print_this_msg("Inside input$brush_select -> Areas", msg_type = "DEBUG")

        brushedPoints(my_val$DT,
                      input$brush_select,
                      xvar = "x",
                      yvar = "y")

        x_min <- input$brush_select$xmin
        x_max <- input$brush_select$xmax
        y_min <- input$brush_select$ymin
        y_max <- input$brush_select$ymax

        if (input$area_unselect) {
          x_dt <- my_val$DT$x[my_val$DT$color == "selected"]
          test_1 <- x_dt > x_min & x_dt < x_max
          y_dt <- my_val$DT$y[my_val$DT$color == "selected"]
          test_2 <- y_dt > y_min & y_dt < y_max
          my_val$DT$color[my_val$DT$color == "selected"][test_1 &
                                                           test_2] <-
            "unselected"
        } else{
          x_dt <- my_val$DT$x[my_val$DT$color == "unselected"]
          test_1 <- x_dt > x_min & x_dt < x_max
          y_dt <- my_val$DT$y[my_val$DT$color == "unselected"]
          test_2 <- y_dt > y_min & y_dt < y_max
          my_val$DT$color[my_val$DT$color == "unselected"][test_1 &
                                                             test_2] <-
            "selected"
        }

        output$main_plot <- shiny::renderPlot({
          draw_plot(input, my_val)
        },
        height = 800,
        width = 800)
      }

    })

    ## Observe event (close_the_path)
    shiny::observeEvent(input$close_the_path, {
      print_this_msg("Inside input$close_the_path event", msg_type = "DEBUG")
      polygons <- my_val$polygon
      if (!is.null(my_val$ggseg)) {
        print_this_msg("Preparing connection.", msg_type = "DEBUG")
        nr_poly <- nrow(polygons)
        polygons[nr_poly,]$xend <- polygons[1,]$x
        polygons[nr_poly,]$yend <- polygons[1,]$y

        print_this_msg("Generating segments", msg_type = "DEBUG")

        ggseg <- ggplot2::geom_segment(
          data = polygons,
          size = 2,
          color = "blue",
          aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend
          )
        )

        # These points cound be inside the polygon:
        x_min <- min(polygons$x)
        x_max <- max(polygons$x)
        y_min <- min(polygons$y)
        y_max <- max(polygons$y)

        my_val$DT$to_be_checked <- 0

        test_1 <-
          my_val$DT$color == "selected" &
          my_val$DT$x > x_min & my_val$DT$x < x_max
        test_2 <-
          my_val$DT$color == "selected" &
          my_val$DT$y > y_min & my_val$DT$y < y_max
        my_val$DT[test_1 & test_2,]$to_be_checked <- 1

        inside_poly <-
          secr::pointsInPolygon(my_val$DT[my_val$DT$to_be_checked == 1, c("x", "y")],
                                polygons[, c("x", "y")])

        my_val$DT[my_val$DT$to_be_checked == 1,]$color[inside_poly] <-
          "unselected"

        my_val$ggseg <- NULL
        my_val$polygon <- NULL
        my_val$DT <-
          my_val$DT[my_val$DT$feature != "POLYGON_COORD",]

        output$main_plot <- shiny::renderPlot({
          draw_plot(input, my_val)
        },
        height = 800,
        width = 800)
      }


    })

    ## Observe event (rem_points_tangents)
    shiny::observeEvent(input$rem_points_tangents, {
      print_this_msg("Removing user points (tangents)", msg_type = "DEBUG")
      my_val$gline <- NULL
      my_val$gpoint <- NULL
      my_val$DT <- my_val$DT[my_val$DT$feature != "POINT_COORDS",]
    })

    ## Observe event (unselect_all)
    shiny::observeEvent(input$unselect_all, {
      print_this_msg("Unselecting all", msg_type = "DEBUG")
      my_val$DT$color <- "unselected"
      output$main_plot <- shiny::renderPlot({
        draw_plot(input, my_val)
      },
      height = 800,
      width = 800)
    })

    ## Observe event (swap_points)
    shiny::observeEvent(input$swap_points, {
      print_this_msg("Swapping points", msg_type = "DEBUG")
      if (input$tool == "Tangents") {
        print_this_msg("Tools: Tangents", msg_type = "DEBUG")
        tmp <- my_val$DT$color == "selected"
        my_val$DT$color <- "selected"
        my_val$DT$color[tmp] <- "unselected"

      } else if (input$tool == "Areas") {
        print_this_msg("Tools: Areass", msg_type = "DEBUG")
        tmp <- my_val$DT$color == "selected"
        my_val$DT$color <- "selected"
        my_val$DT$color[tmp] <- "unselected"

      } else if (input$tool == "Nearest Points") {
        print_this_msg("Tools: Nearest Points", msg_type = "DEBUG")
        tmp <- my_val$DT$color == "selected"
        my_val$DT$color <- "selected"
        my_val$DT$color[tmp] <- "unselected"
      }

      output$main_plot <- shiny::renderPlot({
        draw_plot(input, my_val)
      },
      height = 800,
      width = 800)

    })

    ## Observe event (reset)
    shiny::observeEvent(input$reset, {
      if (input$tool == "Tangents") {
        my_val$DT <- my_val$DT[my_val$DT$feature != "POINT_COORDS",]
        my_val$DT$color <- "selected"
        my_val$gline <- NULL
        my_val$gpoint <- NULL
        my_val$ggseg <- NULL
        my_val$polygon <- NULL


      } else if (input$tool == "Areas") {
        my_val$DT$color <- "selected"

      } else if (input$tool == "Nearest Points") {
        my_val$DT$color <- "selected"
      }

      output$main_plot <- shiny::renderPlot({
        draw_plot(input, my_val)
      },
      height = 800,
      width = 800)
    })

    ## Observe event (download_data)
    output$download_data <- shiny::downloadHandler(
      filename = function() {
        paste("selected_features_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        if (input$save_selected) {
          towrite <- "selected"
        } else{
          towrite <- "unselected"
        }

        write.table(
          my_val$DT[my_val$DT$color == towrite, ],
          file = file,
          sep = "\t",
          col.names = NA,
          quote = FALSE
        )

      }
    )

  }

  shinyApp(
    ui = ui,
    server = server,
    options = list(launch.browser = TRUE)
  )
}

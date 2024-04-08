
# #' Save objects in the current environment
# #'
# #' @param var_to_save
# #' @param path_to_current_ABC_state
# #'
# #' @return
# # ' @export
# #'
# #' @examples
saveEnvir <- function(var_to_save, path_to_current_ABC_state) {
	do.call("save", c(var_to_save, list(file = path_to_current_ABC_state)))
	# save(list=var_to_save, file = path_to_current_ABC_state)
}

# #' Load objects in the current environment of the running ABC
# #'
# #' @param path_to_current_ABC_state
# #'
# #' @return
# # ' @export
# #'
# #' @examples
loadEnvir <- function(path_to_current_ABC_state) {
	load(path_to_current_ABC_state)
}
